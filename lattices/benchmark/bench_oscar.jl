# Oscar side of the Magma/Oscar benchmark: a server that bench.py drives over a pipe.
#
#   julia --project=oscar -t 1 --gcthreads=1 benchmark/bench_oscar.jl CASES
#
# Reads the cases file (bench.py sample), warms up every task on small fixed
# lattices (so that compilation is not timed), prints "@@READY", then answers one
# job per stdin line
#     JOB <case id> <task> <variant> <timeout> <min_time> <max_reps>
# with one line
#     @@RESULT <case id> <task> <variant> <status> <K> <wall,...> <cpu,...> <extra> <value>
# (tab separated; status ok | error | unsupported; times in microseconds per call),
# the same protocol as bench_magma.m.  Anything else on stdout is log output, which
# bench.py ignores.  Julia cannot interrupt a computation, so timeouts and memory
# caps are enforced by bench.py, which kills and restarts this process.
#
# Every timed call gets a fresh lattice (Oscar caches most invariants as attributes
# of the lattice object) and its own Gram matrix U G U^t, the same sequence of bases
# as bench_magma.m uses (see scramble there: Magma reuses results between lattices
# with equal Gram matrices).

const T0 = time()
include(joinpath(@__DIR__, "..", "oscar", "invariants.jl"))
const T_LOAD = time() - T0

struct Case
    id::String
    n::Int
    nplus::Int
    gram::Vector{BigInt}
    gram_iso::Union{Nothing, Vector{BigInt}}
    gram_other::Union{Nothing, Vector{BigInt}}
    theta_prec::Int
end

definite(c::Case) = c.nplus == c.n

# U G U^t for the unimodular U made by max(2n, 12) row operations U[i] += e U[j] drawn from a
# the high bits of a linear congruential generator seeded with k; k = 0 gives G
# itself.  The same function as scramble in bench_magma.m.
function scramble(g::Vector{BigInt}, n::Int, k::Int)
    G = permutedims(reshape(g, n, n))                  # g is row-major
    (k == 0 || n < 2) && return G
    U = BigInt[i == j for i in 1:n, j in 1:n]
    s = k
    for _ in 1:max(2n, 12)
        s = (1103515245 * s + 12345) % 2^31; i = (s ÷ 2^16) % n + 1
        s = (1103515245 * s + 12345) % 2^31; j = (s ÷ 2^16) % (n - 1) + 1
        j >= i && (j += 1)
        s = (1103515245 * s + 12345) % 2^31; e = iseven(s ÷ 2^16) ? 1 : -1
        U[i, :] .+= e .* U[j, :]
    end
    return U * G * permutedims(U)
end

lattice_of(G::Matrix{BigInt}) = integer_lattice(; gram = matrix(QQ, G), cached = false)

# The Gram matrices for the next call of a job, as next_grams in bench_magma.m:
# scramble(gram, k) for the next k whose matrix (and its partner's, for the isometry
# tasks; nothing otherwise) the job has not used yet, or a repeated one after
# MAX_ATTEMPTS tries (a small lattice with a large automorphism group has few).
const MAX_ATTEMPTS = 200
function next_grams(task::String, c::Case, k::Base.RefValue{Int}, used::Set{Matrix{BigInt}})
    partner = task == "isometric_same" ? c.gram_iso : task == "isometric_other" ? c.gram_other : nothing
    for attempt in 1:MAX_ATTEMPTS
        k[] += 1
        G = scramble(c.gram, c.n, k[])
        M = partner === nothing ? nothing : scramble(partner, c.n, k[])
        if attempt == MAX_ATTEMPTS || (!(G in used) && (M === nothing || !(M in used)))
            push!(used, G)
            M === nothing || push!(used, M)
            return G, M
        end
    end
end

parse_gram(s::AbstractString) = s == "-" ? nothing : parse.(BigInt, split(s, ","))

function read_cases(path::String)
    cases = Dict{String, Case}()
    lines = readlines(path)
    header = split(lines[1], "\t")
    col(name) = findfirst(==(name), header)
    for line in lines[2:end]
        isempty(line) && continue
        f = split(line, "\t")
        c = Case(String(f[col("id")]), parse(Int, f[col("rank")]), parse(Int, f[col("nplus")]),
                 parse_gram(f[col("gram")]), parse_gram(f[col("gram_iso")]),
                 parse_gram(f[col("gram_other")]), parse(Int, f[col("theta_prec")]))
        cases[c.id] = c
    end
    return cases
end

# ---------------------------------------------------------------------------
# Tasks.  Each implementation is (applies, setup, run): applies(case) says whether
# it can run at all (otherwise the job is "unsupported"); setup(L, case, M) does the
# untimed preparation on the fresh lattice L (M: the partner's Gram matrix for the
# isometry tasks, see next_grams) and returns a state;
# run(L, state, case) is timed and returns the value as a string in the format the
# Magma server uses, or nothing for "unsupported" (decided only at run time, e.g. by
# the group order).  Variant "default" is the one compared against Magma's default.

struct Impl
    applies::Function
    setup::Function
    run::Function
end

always(c) = true
nosetup(L, c, M) = nothing

const TASKS = Dict{String, Dict{String, Impl}}(
    "class_number" => Dict(
        # genus_representatives dispatches like the pipeline's GenusReps: rank 2 via
        # binary forms, definite rank >= 3 by Kneser neighbours (mass-checked),
        # indefinite rank >= 3 by spinor genera.
        "default" => Impl(always, nosetup, (L, s, c) -> string(length(genus_representatives(L)))),
        "orbit" => Impl(c -> definite(c) && c.n >= 3, nosetup,
                        (L, s, c) -> string(length(enumerate_definite_genus(L, :orbit)))),
        "random" => Impl(c -> definite(c) && c.n >= 3, nosetup,
                         (L, s, c) -> string(length(enumerate_definite_genus(L, :random)))),
    ),
    "mass" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> qq_string(mass(L))),
    ),
    "aut_order" => Dict(
        # isometry_group is what the pipeline needs (aut_group, aut_label, is_chiral);
        # automorphism_group_order skips building the GAP matrix group.
        "default" => Impl(definite, nosetup, (L, s, c) -> string(order(isometry_group(L)))),
        "order_only" => Impl(definite, nosetup, (L, s, c) -> string(automorphism_group_order(L))),
    ),
    "aut_label" => Dict(
        "default" => Impl(definite, (L, c, M) -> isometry_group(L), (L, G, c) -> aut_label(G)),
    ),
    "minimum" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> qq_string(minimum(L))),
    ),
    "kissing" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> string(kissing_number(L))),
    ),
    "successive_minima" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> successive_minima_string(L)),
    ),
    "root_system" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> root_string(L)),
    ),
    "theta" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> join(string.(theta_series(L, c.theta_prec)), ",")),
    ),
    "orth_decomp" => Dict(
        "default" => Impl(definite, nosetup, (L, s, c) -> decomposition_string(irreducible_components(L))),
    ),
    "isometric_same" => Dict(
        "default" => Impl(c -> c.gram_iso !== nothing, (L, c, M) -> lattice_of(M),
                          (L, M, c) -> string(is_isometric(L, M))),
    ),
    "isometric_other" => Dict(
        "default" => Impl(c -> c.gram_other !== nothing, (L, c, M) -> lattice_of(M),
                          (L, M, c) -> string(is_isometric(L, M))),
    ),
    # |O(q_L)| alone (the numerator of the Festi-Veniani index); Magma has no
    # orthogonal groups of torsion quadratic forms.
    "oq_order" => Dict(
        "default" => Impl(always, nosetup, (L, s, c) -> string(order(orthogonal_group(discriminant_group(L))))),
    ),
    "disc_image" => Dict(
        # The isometry group is computed in the setup, as the Magma side computes the
        # automorphism group there; the timed part is the image in O(q_L).
        "default" => Impl(always, (L, c, M) -> (definite(c) && isometry_group(L); nothing),
                          (L, s, c) -> begin
                              im, oq = discriminant_image(L)
                              EXTRA[] = "Oq=$(oq)"
                              string(im)
                          end),
    ),
)

const EXTRA = Ref("")       # task-specific extra output (e.g. |O(q)|), reset per job

cpu_seconds() = ccall(:clock, Clong, ()) / 1_000_000    # CLOCKS_PER_SEC on Linux/glibc

clean(s) = replace(string(s), r"[\t\r\n]+" => " ")

# Timing as in bench_magma.m, so the two are measured the same way: blocks of K
# calls, each on a fresh lattice with its untimed setup done beforehand, K growing
# tenfold until a block takes BLOCK seconds (or preparing one takes PREP_LIMIT);
# then blocks of the final K until min_time seconds or max_reps blocks.  Julia's
# clock is finer than Magma's 10 ms, but GC pauses make single fast calls noisy too.
const BLOCK = 0.2
const PREP_LIMIT = 5.0
kmax(n) = n <= 2 ? 1000 : 10_000     # rank 2 has few small bases (see next_grams)

struct Block
    values::Vector{Any}
    wall::Float64
    cpu::Float64
    prep::Float64
    gc::Float64
    compile::Float64
    bytes::Int
end

# K lattices on fresh Gram matrices (next_grams), with their setups
function block(task::String, c::Case, impl::Impl, K::Int, k, used)
    t = time()
    Ls = ZZLat[]
    sts = Any[]
    for _ in 1:K
        G, M = next_grams(task, c, k, used)
        L = lattice_of(G)
        push!(Ls, L)
        push!(sts, impl.setup(L, c, M))
    end
    prep = time() - t
    GC.gc()        # untimed: a block pays for its own garbage, not for earlier jobs'
    c0 = cpu_seconds()
    s = @timed [impl.run(Ls[i], sts[i], c) for i in 1:K]
    return Block(s.value, s.time, cpu_seconds() - c0, prep, s.gctime, s.compile_time, s.bytes)
end

const UNSUPPORTED = ("unsupported", 0, Float64[], Float64[], "", "")

function run_job(task::String, c::Case, impl::Impl; min_time = 0.5, max_reps = 5)
    impl.applies(c) || return UNSUPPORTED
    EXTRA[] = ""
    K = 1
    k = Ref(-1)                   # the first call gets k = 0, the Gram matrix as given
    used = Set{Matrix{BigInt}}()
    b = block(task, c, impl, K, k, used)
    # A call that compiled code paths the warm-up missed: time another (up to twice).
    for _ in 1:2
        b.compile > 0.1 * b.wall || break
        b = block(task, c, impl, K, k, used)
    end
    value = b.values[1]
    value === nothing && return UNSUPPORTED
    check(b) = all(==(value), b.values) || error("value changed between calls: $value, then $(first(filter(!=(value), b.values)))")
    while b.wall < BLOCK && K < kmax(c.n) && b.prep < PREP_LIMIT
        K = min(kmax(c.n), 10K)
        b = block(task, c, impl, K, k, used)
        check(b)
    end
    blocks = [b]
    while length(blocks) < max_reps && sum(x.wall for x in blocks) < min_time
        b = block(task, c, impl, K, k, used)
        check(b)
        push!(blocks, b)
    end
    calls = K * length(blocks)
    extra = join(filter(!isempty, [EXTRA[],
                 "gc=$(round(Int, 1e6 * sum(x.gc for x in blocks) / calls))",
                 "compile=$(round(Int, 1e6 * sum(x.compile for x in blocks) / calls))",
                 "alloc=$(sum(x.bytes for x in blocks) ÷ calls)"]), ";")
    return ("ok", K, [x.wall / K for x in blocks], [x.cpu / K for x in blocks], extra, value)
end

usec(ts) = isempty(ts) ? "-" : join([string(round(Int, 1e6 * t)) for t in ts], ",")

function reply(id, task, variant, status, K, walls, cpus, extra, value)
    println(join(["@@RESULT", id, task, variant, status, string(K), usec(walls), usec(cpus),
                  isempty(extra) ? "-" : clean(extra), isempty(value) ? "-" : clean(value)], "\t"))
    flush(stdout)
end

# Small fixed lattices covering the code paths (definite rank 2, 3, 5 and 8, odd and
# even, class numbers 1 and 3; indefinite rank 2 and 3), run through every task
# before the first job.  The same ones as in bench_magma.m.
function warm_cases()
    g(v) = BigInt.(v)
    A2_6 = [2, -1, 0, -1, 2, 0, 0, 0, 6]
    GE8 = gram_matrix(root_lattice(:E, 8))
    E8 = [Int(ZZ(GE8[i, j])) for i in 1:8 for j in 1:8]
    return [
        Case("warm1", 3, 3, g(A2_6), g([2, 1, 0, 1, 2, 0, 0, 0, 6]), g([2, 0, 0, 0, 2, 0, 0, 0, 6]), 20),
        Case("warm2", 8, 8, g(E8), g(E8), nothing, 6),
        Case("warm3", 2, 2, g([1, 0, 0, 5]), g([1, 1, 1, 6]), g([2, 1, 1, 3]), 20),
        Case("warm4", 3, 3, g([1, 0, 0, 0, 1, 0, 0, 0, 3]), nothing, nothing, 20),
        Case("warm5", 2, 1, g([2, 1, 1, -2]), g([2, -1, -1, -2]), nothing, 0),
        Case("warm6", 3, 2, g([2, 1, 0, 1, 2, 0, 0, 0, -18]), nothing, nothing, 0),
        Case("warm7", 3, 3, g([2, 0, 0, 0, 2, -1, 0, -1, 16]), g([2, 0, 0, 0, 2, 1, 0, 1, 16]), nothing, 40),
        Case("warm8", 5, 5, g([2, -1, 1, -1, 1, -1, 2, -1, 0, -1, 1, -1, 4, 1, 2, -1, 0, 1, 4, -1, 1, -1, 2, -1, 6]),
             nothing, nothing, 12),
    ]
end

function main(path::String)
    cases = read_cases(path)
    t1 = time()
    for c in warm_cases(), (task, impls) in TASKS, (variant, impl) in impls
        try
            if impl.applies(c)                  # one untimed call, to compile it
                L = lattice_of(scramble(c.gram, c.n, 0))
                M = task == "isometric_same" ? c.gram_iso : c.gram_other
                impl.run(L, impl.setup(L, c, M === nothing ? nothing : scramble(M, c.n, 0)), c)
            end
        catch e
            println(stderr, "warm-up $(c.id) $task/$variant: ", first(clean(sprint(showerror, e)), 300))
        end
    end
    t_warm = time() - t1
    versions = "Oscar $(Oscar.VERSION_NUMBER) Hecke $(pkgversion(Oscar.Hecke)) Julia $(VERSION)"
    println()                 # end any progress line the warm-up left unterminated
    println("@@READY\t$(round(T_LOAD; digits = 1))\t$(round(t_warm; digits = 1))\t$versions")
    flush(stdout)
    for line in eachline(stdin)
        f = split(line, "\t")
        # JOB <case id> <task> <variant> <timeout> <min_time> <max_reps>; the timeout
        # is enforced by bench.py
        (length(f) >= 7 && f[1] == "JOB") || continue
        id, task, variant = String(f[2]), String(f[3]), String(f[4])
        min_time = parse(Float64, f[6]); max_reps = parse(Int, f[7])
        c = get(cases, id, nothing)
        impl = get(get(TASKS, task, Dict{String, Impl}()), variant, nothing)
        if c === nothing || impl === nothing
            reply(id, task, variant, "error", 0, Float64[], Float64[], "", "unknown case, task or variant")
            continue
        end
        try
            reply(id, task, variant, run_job(task, c, impl; min_time, max_reps)...)
        catch e
            reply(id, task, variant, "error", 0, Float64[], Float64[], "", first(clean(sprint(showerror, e)), 300))
        end
    end
end

main(ARGS[1])
