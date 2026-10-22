# Resolve the label lookups connect could not: lattices in an indefinite genus of
# rank >= 3 with several classes.  Those classes are the genus's spinor genera, and
# Magma has no way to tell which one a given lattice is in (FindLabel in
# connect_genus.m returns null and the genus), but Oscar's is_isometric handles
# indefinite lattices.  Each line `column|genus label|gram` of
# indefinite_lookups/<rank>/<nplus>/<source label> is a lattice whose label belongs
# in that column of the source lattice's advanced record; this finds the class of the
# genus isometric to it and writes the class's label there.
#
# Usage, from lattices/ after connect:
#   julia --project=oscar oscar/resolve_indefinite.jl [root]
# Exits 1 if some lattice matches no class or several (an enumeration error); a genus
# whose classes are not all stored leaves the column null.

include(joinpath(@__DIR__, "load_local.jl"))

function resolve(root::String)
    fmts = pipeline_formats(root)
    adv = read_format(joinpath(root, "lat_advanced.format"))
    dir = joinpath(root, "indefinite_lookups")
    if !isdir(dir)
        println("no indefinite lookups")
        return 0
    end
    classes = Dict{String, Union{Nothing, Vector{ZZLat}}}()
    resolved = skipped = 0
    problems = String[]
    for (d, _, files) in walkdir(dir), src in files
        path = label_path(root, "lattice_advanced_data", src)
        if !isfile(path)
            push!(problems, "$src: no advanced record")
            continue
        end
        lines = split(read(path, String), "\n")
        fields = String.(split(lines[1], "|"))            # split keeps empty fields
        length(fields) == length(adv) || error("$path: $(length(fields)) fields, lat_advanced.format has $(length(adv))")
        changed = false
        for line in eachline(joinpath(d, src))
            isempty(line) && continue
            col, genus, gram = split(line, "|")
            ci = findfirst(==(col), adv)
            ci === nothing && error("$src: unknown column $col")
            Ls = get!(classes, genus) do
                meta = genus_meta(root, fmts, genus)
                meta === nothing ? nothing : load_class_lattices(root, fmts, meta)
            end
            if Ls === nothing                              # not all classes stored
                skipped += 1
                continue
            end
            entries = parse.(BigInt, split(strip(gram, ['{', '}']), ","))
            n = isqrt(length(entries))
            Q = integer_lattice(gram = matrix(ZZ, n, n, entries))
            hits = [L for L in Ls if is_isometric(Q, L)]
            if length(hits) != 1
                push!(problems, "$src $col: $(length(hits)) classes of $genus are isometric to it")
                continue
            end
            fields[ci] = get_attribute(hits[1], :lmfdb_label)
            changed = true
            resolved += 1
        end
        if changed
            lines[1] = join(fields, "|")
            tmp = path * ".resolve.tmp"
            write(tmp, join(lines, "\n"))
            mv(tmp, path; force = true)
        end
    end
    println("indefinite lookups: $resolved resolved, $skipped left null (genus not fully stored), ",
            "$(length(problems)) problems")
    foreach(p -> println("  ", p), problems)
    return isempty(problems) ? 0 : 1
end

exit(resolve(length(ARGS) >= 1 ? ARGS[1] : "."))
