# A pool of persistent Oscar workers (Distributed.jl) that runs a per-genus task
# over a list of labels.
#
# Why a pool: a fresh Julia process pays ~30 s to load Oscar (more under load), and
# the first call of each kind of computation used to pay minutes of compilation
# (110 s for the first definite `representatives` at load ~160).  The package's
# precompile workload (workload.jl) caches that compiled code on disk, and the pool
# pays the load once per worker instead of once per job.
#
# Why kill-and-replace: Julia cannot interrupt a computation, so a task that
# overruns its time limit or memory cap is stopped by killing its worker process,
# which is then replaced; the other workers keep going.
#
# Work goes out in batches (one remote call per batch, so cheap genera do not make
# the master a bottleneck).  The worker writes each label's record itself, atomically,
# as it finishes it, and works through the batch in order, so the master can see which
# label a stuck worker is on: the first one without a record.  That label alone is
# charged the timeout (a failure record), and the rest of its batch goes back on the
# queue.

using Distributed

function out_path(out_root::AbstractString, label::AbstractString)
    parts = split(label, ".")
    return joinpath(out_root, String(parts[1]), String(parts[2]), String(label))
end

record_status(path) = isfile(path) ? (m = match(r"^status=(\S+)", first(eachline(path))); m === nothing ? "bad" : m[1]) : nothing

function write_record(path::AbstractString, text::AbstractString)
    mkpath(dirname(path))
    tmp = string(path, ".tmp.", getpid())
    write(tmp, text)
    mv(tmp, path; force = true)
end

# Labels "selftest.<kind>.<n>" exercise the pool itself (oscar/test_pool.jl): sleep
# n seconds (yielding, or in C without yielding), allocate n GB and hold it, or throw.
function selftest_task(label::AbstractString, root::AbstractString)
    _, kind, n = split(label, ".")
    if kind == "sleep"
        sleep(parse(Float64, n))
    elseif kind == "csleep"
        # blocks in C without yielding, like a long Oscar/GAP call: the worker cannot
        # answer the master until it returns
        ccall(:sleep, Cuint, (Cuint,), parse(Int, n))
    elseif kind == "alloc"
        x = ones(UInt8, round(Int, parse(Float64, n) * 2^30))
        sleep(30)
        return "status=ok seconds=0 bytes=$(sum(x))\n"
    elseif kind == "error"
        error("selftest error $n")
    end
    return "status=ok seconds=0\n"
end

const TASKS = Dict{String, Function}("enumerate" => enumerate_task, "resolve" => resolve_task, "embed" => embed_task,
                                     "selftest" => selftest_task)

"""
    run_batch(task, labels, root, out_root) -> Int

Worker side: run `task` on each label in order, writing each record as soon as it is
done (an error becomes an error record).  Returns how many records were written.
"""
function run_batch(task::String, labels::Vector{String}, root::String, out_root::String)
    f = TASKS[task]
    n = 0
    for label in labels
        p = out_path(out_root, label)
        record_status(p) == "ok" && continue
        t0 = time()
        text = try
            f(label, root)
        catch e
            failure_record("error", time() - t0; message = sprint(showerror, e))
        end
        write_record(p, text)
        n += 1
    end
    return n
end

function rss_bytes(pid::Integer)
    statm = "/proc/$pid/statm"
    if isfile(statm)
        pages = parse(Int, split(read(statm, String))[2])
        return pages * 4096
    end
    out = try read(`ps -o rss= -p $pid`, String) catch; "" end      # macOS: no /proc
    return isempty(strip(out)) ? 0 : 1024 * parse(Int, strip(out))
end

mutable struct PoolWorker
    id::Int
    pid::Int
end

logp(msg...) = (println(stdout, "pool: ", msg...); flush(stdout))

function start_worker(exeflags)
    # A worker that has not reported back within JULIA_WORKER_TIMEOUT (default 60 s)
    # makes addprocs throw; under load ~200 workers took ~100 s to start.
    haskey(ENV, "JULIA_WORKER_TIMEOUT") || (ENV["JULIA_WORKER_TIMEOUT"] = "900")
    for attempt in 1:3
        try
            w = only(addprocs(1; exeflags = exeflags))
            remotecall_eval(Main, [w], :(using LMFDBLattices))
            return PoolWorker(w, remotecall_fetch(getpid, w))
        catch e
            logp("starting a worker failed (attempt $attempt): ", sprint(showerror, e))
            attempt == 3 && rethrow()
        end
    end
end

function kill_worker(pw::PoolWorker)
    # Mark the worker as being removed before killing it: otherwise Distributed reads
    # the dead socket as an unexpected death and prints "Worker N terminated." and an
    # unhandled EOFError.  rmprocs would mark it too, but only once it gets
    # worker_lock, which addprocs holds while a replacement worker starts.
    try Distributed.set_worker_state(Distributed.worker_from_id(pw.id), Distributed.W_TERMINATING) catch end
    try run(`kill -9 $(pw.pid)`) catch end
    try rmprocs(pw.id; waitfor = 0) catch end
end

"""
    run_pool(task, labels; root, out_root, jobs, timeout, memlimit_gb, batch, retry_failed, ledger)

Run `task` over `labels` on `jobs` workers.  Labels with an `ok` record are skipped
(resume); with `retry_failed`, so are only those, otherwise any existing record is
kept.  `timeout` is per label (seconds), `memlimit_gb` per worker process.  Appends
one line per batch to `ledger` (tab-separated: worker, labels, done, seconds, outcome).
Returns a Dict of record statuses.
"""
function run_pool(task::String, labels::Vector{String}; root::String = ".", out_root::String = "oscar_genus",
                  jobs::Int = 4, timeout::Real = 600, memlimit_gb::Real = 16, batch::Int = 20,
                  retry_failed::Bool = false, ledger::String = "", project::String = dirname(Base.active_project()))
    haskey(TASKS, task) || error("unknown task $task (known: $(join(keys(TASKS), ", ")))")
    todo = String[]
    for l in labels
        s = record_status(out_path(out_root, l))
        (s === nothing || (retry_failed && s != "ok")) && push!(todo, l)
    end
    logp("$(length(labels)) labels, $(length(labels) - length(todo)) already done, $(length(todo)) to run on $jobs workers")
    isempty(todo) && return summarize(labels, out_root)
    queue = Channel{Vector{String}}(Inf)
    for k in 1:batch:length(todo)
        put!(queue, todo[k:min(k + batch - 1, end)])
    end
    hint = max(1, floor(Int, 0.75 * memlimit_gb))
    exeflags = `--project=$project -t 1 --heap-size-hint=$(hint)G --startup-file=no`
    memlimit = memlimit_gb * 2^30
    pending = Ref(length(todo))                 # labels not yet finished or failed
    lk = ReentrantLock()
    log_line(fields...) = isempty(ledger) || lock(lk) do
        open(ledger, "a") do io; println(io, join(fields, "\t")); end
    end
    t_start = time()
    starting = [@async start_worker(exeflags) for _ in 1:jobs]
    pws = fetch.(starting)
    logp("$(length(pws)) workers up after $(round(time() - t_start; digits = 1)) s")
    function drive(pw::PoolWorker)
        while pending[] > 0
            b = try take!(queue) catch; break end
            try
                drive_batch(pw, b)
            catch e
                # Anything unexpected here must not kill this worker's loop silently:
                # its batch would never be counted and the pool would wait forever.
                logp("worker $(pw.id) failed on a batch of $(length(b)): ", sprint(showerror, e, catch_backtrace()))
                kill_worker(pw)
                for l in b
                    record_status(out_path(out_root, l)) === nothing &&
                        write_record(out_path(out_root, l), failure_record("crashed", 0.0; message = "pool error"))
                end
                lock(lk) do; pending[] -= length(b); end
                newpw = start_worker(exeflags)
                pw.id, pw.pid = newpw.id, newpw.pid
            end
            pending[] <= 0 && close(queue)
        end
    end
    function drive_batch(pw::PoolWorker, b::Vector{String})
            # remotecall_fetch in a local task, polled with istaskdone: isready on a
            # remotecall Future asks the worker, and a worker inside a long Oscar/GAP
            # call does not answer until it returns -- so the timeout check never ran.
            # (the task returns an exception rather than throwing it: after a kill its
            # fetch fails, and nobody waits on it any more)
            fut = @async try remotecall_fetch(run_batch, pw.id, task, b, root, out_root) catch e; e end
            k, t_label, t0 = 1, time(), time()
            outcome = "ok"
            while !istaskdone(fut)
                sleep(0.5)
                while k <= length(b) && record_status(out_path(out_root, b[k])) !== nothing
                    k += 1; t_label = time()
                end
                if time() - t_label > timeout
                    outcome = "timeout"
                elseif rss_bytes(pw.pid) > memlimit
                    outcome = "memory"
                end
                outcome == "ok" || break
            end
            if outcome == "ok"
                e = fetch(fut)
                if e isa Exception
                    outcome = e isa ProcessExitedException || e isa EOFError ||
                              (e isa RemoteException && e.captured.ex isa ProcessExitedException) ? "crashed" : "error"
                end
            end
            while k <= length(b) && record_status(out_path(out_root, b[k])) !== nothing
                k += 1
            end
            done = k - 1
            if outcome == "ok"
                lock(lk) do; pending[] -= length(b); end
            else
                # the label being worked on is charged with the failure; the rest of
                # the batch goes back on the queue
                outcome == "error" || kill_worker(pw)
                stuck = b[min(k, end)]
                if record_status(out_path(out_root, stuck)) === nothing
                    write_record(out_path(out_root, stuck), failure_record(outcome, time() - t_label))
                end
                rest = b[k + 1:end]
                isempty(rest) || put!(queue, rest)
                lock(lk) do; pending[] -= length(b) - length(rest); end
                logp("worker $(pw.id): $(stuck) $(outcome) after $(round(time() - t_label; digits = 1)) s")
                if outcome != "error"
                    newpw = start_worker(exeflags)
                    pw.id, pw.pid = newpw.id, newpw.pid
                end
            end
            log_line(pw.id, length(b), done, round(time() - t0; digits = 2), outcome)
    end
    @sync for pw in pws
        @async drive(pw)
    end
    for pw in pws
        try rmprocs(pw.id; waitfor = 30) catch end
        kill_worker(pw)                       # rmprocs does not stop a busy worker
    end
    logp("finished in $(round(time() - t_start; digits = 1)) s")
    return summarize(labels, out_root)
end

function summarize(labels, out_root)
    counts = Dict{String, Int}()
    for l in labels
        s = something(record_status(out_path(out_root, l)), "missing")
        counts[s] = get(counts, s, 0) + 1
    end
    logp(join(["$k=$v" for (k, v) in sort(collect(counts))], " "))
    return counts
end

function parse_cli(args)
    opts = Dict{String, String}()
    flags = Set(["retry-failed"])
    i = 1
    while i <= length(args)
        startswith(args[i], "--") || error("unexpected argument $(args[i])")
        key = args[i][3:end]
        if key in flags
            opts[key] = "true"; i += 1
        else
            i + 1 <= length(args) || error("missing value for --$key")
            opts[key] = args[i + 1]; i += 2
        end
    end
    return opts
end

"""
    main_pool(ARGS)

Command line (oscar/pool.jl):
  --task enumerate --labels-file F [--jobs 4] [--timeout 600] [--memlimit 16]
  [--batch 20] [--data-root .] [--out-root oscar_genus] [--ledger FILE] [--retry-failed]
Exits 1 if any label ended without an ok record.
"""
function main_pool(args)
    o = parse_cli(args)
    labels = filter(!isempty, strip.(readlines(o["labels-file"])))
    counts = run_pool(o["task"], String.(labels);
                      root = get(o, "data-root", "."), out_root = get(o, "out-root", "oscar_genus"),
                      jobs = parse(Int, get(o, "jobs", "4")), timeout = parse(Float64, get(o, "timeout", "600")),
                      memlimit_gb = parse(Float64, get(o, "memlimit", "16")), batch = parse(Int, get(o, "batch", "20")),
                      retry_failed = haskey(o, "retry-failed"), ledger = get(o, "ledger", ""))
    return all(k == "ok" for k in keys(counts)) ? 0 : 1
end
