# Assemble lat_inclusions from the raw per-lattice files that driver.jl wrote
# under <inclusions-root>/data/<jobtag>/<rank>/<nplus>/<domain_label>/.  Port of Stevell's
# writing_file.jl with three fixes: ranks are parsed from labels rather than
# from the last character of "rankN" (which breaks at rank 10), the K3 lattice
# label is our scheme's 22.19.1.0.1 rather than the stale hard-coded one, and
# the output is real psycodict copy format (three header lines, then exactly
# one row per line -- no blank line between rows).
#
# Usage: julia oscar/write_table.jl --inclusions-root inclusions --out lat_inclusions.tbl
#        [--k3-label 22.19.1.0.1] [--lattice-root lattice_basic_data]
# Plain Julia, no Oscar needed.  With --lattice-root, every referenced label is
# checked against the generated lattice files and dangling references are
# reported (they make the row unloadable, so they fail the run).

const LAT_INCLUSIONS_HEADER = (
  ("domain", "text"),
  ("codomain", "text"),
  ("primitive", "boolean"),
  ("quotient_invs", "integer[]"),
  ("domain_rank", "smallint"),
  ("codomain_rank", "smallint"),
  ("quotient_rank", "smallint"),
  ("multiplicity", "integer"),
  ("orthogonal_complement", "numeric[]"),
  ("quotient_gens", "numeric[]"),
  ("images", "numeric[]"),
)

label_rank(label) = parse(Int, first(split(label, ".")))

# "{a,b},{c,d}" nesting: an all-empty list collapses to {} (postgres arrays
# cannot mix empty and nonempty members, and indefinite codomains carry no
# complement coordinates)
function nest(vecs::Vector{String})
    all(v -> v == "{}", vecs) && return "{}"
    return "{" * join(vecs, ",") * "}"
end

braces(s) = "{" * strip(s, ['[', ']']) * "}"

function parse_args(args)
    opts = Dict{String, String}("inclusions-root" => "inclusions",
                                "out" => "lat_inclusions.tbl",
                                "k3-label" => "22.19.1.0.1")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") || error("unexpected argument $(args[i])")
        key = args[i][3:end]
        opts[key] = args[i + 1]
        i += 2
    end
    return opts
end

function main()
    opts = parse_args(ARGS)
    root = opts["inclusions-root"]
    k3_label = opts["k3-label"]
    latroot = get(opts, "lattice-root", "")
    rows = 0
    dangling = String[]

    function check_label(label)
        isempty(latroot) && return
        parts = split(label, ".")
        p = joinpath(latroot, String(parts[1]), String(parts[2]), String(label))
        isfile(p) || push!(dangling, String(label))
    end

    open(opts["out"], "w") do out
        write(out, join(first.(LAT_INCLUSIONS_HEADER), "|") * "\n")
        write(out, join(last.(LAT_INCLUSIONS_HEADER), "|") * "\n\n")
        # Only committed jobs are read: the driver stages under staging/<jobtag>/
        # and renames to data/<jobtag>/ (manifest included) atomically on clean
        # completion, so partial or duplicated output from killed/retried jobs
        # never reaches the table.
        data = joinpath(root, "data")
        isdir(data) || return
        for jdir in sort(readdir(data; join = true)), rdir in
                (isdir(jdir) ? sort(readdir(jdir; join = true)) : String[])
            (isdir(rdir) && all(isdigit, basename(rdir))) || continue
            for sdir in sort(readdir(rdir; join = true))
            isdir(sdir) || continue
            for ddir in sort(readdir(sdir; join = true))
                isdir(ddir) || continue
                d_label = basename(ddir)
                r = label_rank(d_label)
                check_label(d_label)
                for fname in sort(readdir(ddir))
                    path = joinpath(ddir, fname)
                    if fname == "k3_embed.txt"
                        count = strip(read(path, String))
                        check_label(k3_label)
                        write(out, join([d_label, k3_label, "t", "{}", "$r", "22",
                                         "$(22 - r)", count, "{}", "{}", "{}"], "|") * "\n")
                        rows += 1
                    elseif fname == "overlat.txt"
                        groups = Dict{Tuple{String, String}, Vector{String}}()
                        for line in eachline(path)
                            c_label, qi, qg = split(line, "_")
                            push!(get!(groups, (String(c_label), "{$qi}"), String[]),
                                  braces(qg))
                        end
                        for ((c_label, qi), qgs) in sort(collect(groups))
                            check_label(c_label)
                            write(out, join([d_label, c_label, "f", qi, "$r",
                                             "$(label_rank(c_label))", "0",
                                             "$(length(qgs))", "{}", nest(qgs), "{}"],
                                            "|") * "\n")
                            rows += 1
                        end
                    elseif fname == "prim_emb.txt"
                        groups = Dict{String, Vector{String}}()
                        for line in eachline(path)
                            c_label, oc = split(line, "_")
                            push!(get!(groups, String(c_label), String[]), braces(oc))
                        end
                        for (c_label, ocs) in sort(collect(groups))
                            check_label(c_label)
                            write(out, join([d_label, c_label, "t", "{}", "$r",
                                             "$(label_rank(c_label))",
                                             "$(label_rank(c_label) - r)",
                                             "$(length(ocs))", nest(ocs), "{}", "{}"],
                                            "|") * "\n")
                            rows += 1
                        end
                    end
                end
            end
            end
        end
    end
    println("wrote $rows rows to $(opts["out"])")
    if !isempty(dangling)
        println("ERROR: $(length(dangling)) referenced labels have no lattice file, e.g. ",
                join(unique(dangling)[1:min(5, end)], ", "))
        exit(1)
    end
end

main()
