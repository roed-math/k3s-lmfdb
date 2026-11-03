# Fill in festi_veniani_index = [O(q_L) : image of O(L)], the Festi-Veniani index, in
# lattice_basic_data.  O(q_L) is the orthogonal group of the discriminant form; Magma has
# none for torsion quadratic forms, so fill computes only the order of the image of O(L)
# on the discriminant group (which lies in O(q_L)) and records it in
#     festi_veniani/<rank>/<nplus>/<genus label>:  "<lattice label>|<image order>" lines,
# leaving the column \N.  This computes |O(q_L)| once per genus with Oscar (GAP), from the
# genus representative in genera_basic, and writes |O(q_L)| / |image| into each lattice's
# record.  The sidecar holds the image order, so a rerun gives the same values.
#
# --legacy also converts data filled before 2026-11 (no sidecar), whose column holds
# #Aut(L^#/L) / |image| with Aut(L^#/L) the automorphism group of the abelian group: it
# recovers |image| = #Aut(L^#/L) / stored value, writes the sidecar first, then patches.
#
# Usage, from lattices/:
#   julia --project=oscar oscar/fix_festi_veniani.jl [--legacy] ROOT [GENUS_LABELS_FILE ...]
# Without label files every genus with a sidecar (and, with --legacy, every genus with
# stored lattices) is processed.  Keep festi_veniani/ with the data: it is what makes
# the stored values interpretable.  Exits 1 if some genus is inconsistent
# (an image order that does not divide |O(q_L)|, or a legacy value that does not divide
# #Aut(L^#/L)); those genera are left as they were.

include(joinpath(@__DIR__, "load_local.jl"))

"""
    abelian_aut_order(invs) -> ZZRingElem

The order of the automorphism group of the finite abelian group with invariant factors
`invs`, by the formula of Hillar and Rhea ("Automorphisms of finite abelian groups",
Amer. Math. Monthly 114 (2007)) on each p-part: for exponents e_1 <= ... <= e_n,
prod_k (p^{d_k} - p^{k-1}) * prod_j p^{e_j (n - d_j)} * prod_i p^{(e_i - 1)(n - c_i + 1)}
with d_k = max{l : e_l = e_k}, c_k = min{l : e_l = e_k}.
"""
function abelian_aut_order(invs::Vector{ZZRingElem})
    exps = Dict{ZZRingElem, Vector{Int}}()
    for m in invs, (p, e) in factor(m)
        push!(get!(exps, p, Int[]), e)
    end
    total = ZZ(1)
    for (p, es) in exps
        sort!(es)
        n = length(es)
        d = [findlast(==(es[k]), es) for k in 1:n]
        c = [findfirst(==(es[k]), es) for k in 1:n]
        for k in 1:n
            total *= p^d[k] - p^(k - 1)
            total *= p^(es[k] * (n - d[k]))
            total *= p^((es[k] - 1) * (n - c[k] + 1))
        end
    end
    return total
end

parse_list(s::AbstractString) = [ZZ(parse(BigInt, x)) for x in split(strip(s, ['{', '}', ' ']), ",") if !isempty(strip(x))]

function genus_fields(root::String, fmt::Vector{String}, label::AbstractString)
    path = label_path(root, "genera_basic", label)
    isfile(path) || return nothing
    pieces = split(first(eachline(path)), "|")
    length(pieces) == length(fmt) || return nothing
    return Dict(zip(fmt, pieces))
end

function genus_lattice(g::Dict)
    entries = parse_list(g["rep"])
    n = isqrt(length(entries))
    return integer_lattice(; gram = matrix(QQ, n, n, entries))
end

read_sidecar(path::String) =
    [(String(a), ZZ(parse(BigInt, b))) for (a, b) in (split(l, "|") for l in eachline(path) if !isempty(l))]

function write_atomically(path::String, text::String)
    tmp = path * ".fv.tmp"
    write(tmp, text)
    mv(tmp, path; force = true)
end

# The lattice_basic_data records of a genus: label.1, label.2, ... while they exist.
function genus_records(root::String, label::AbstractString)
    out = String[]
    i = 1
    while isfile(label_path(root, "lattice_basic_data", "$label.$i"))
        push!(out, "$label.$i")
        i += 1
    end
    return out
end

# Legacy data: recover the image orders from the stored #Aut(L^#/L) / |image|.  Returns
# the sidecar entries (also written to disk), [] if there is nothing to convert, or a
# problem string.
function legacy_sidecar(root, lfmt, fvi, label, g, sidecar)
    autA = abelian_aut_order(parse_list(g["discriminant_group_invs"]))
    entries = Tuple{String, ZZRingElem}[]
    for lab in genus_records(root, label)
        v = split(first(eachline(label_path(root, "lattice_basic_data", lab))), "|")[fvi]
        v in ("\\N", "None", "") && continue
        stored = ZZ(parse(BigInt, v))
        (stored > 0 && divides(autA, stored)[1]) ||
            return "$lab: legacy value $stored does not divide #Aut(L^#/L) = $autA"
        push!(entries, (lab, div(autA, stored)))
    end
    if !isempty(entries)
        mkpath(dirname(sidecar))
        write_atomically(sidecar, join(["$a|$b" for (a, b) in entries], "\n") * "\n")
    end
    return entries
end

function process_genus(root, gfmt, lfmt, fvi, label; legacy::Bool)
    sidecar = label_path(root, "festi_veniani", label)
    g = genus_fields(root, gfmt, label)
    g === nothing && return (0, "$label: no genera_basic record")
    if isfile(sidecar)
        entries = read_sidecar(sidecar)
    elseif legacy
        entries = legacy_sidecar(root, lfmt, fvi, label, g, sidecar)
        entries isa String && return (0, entries)
    else
        return (0, nothing)
    end
    isempty(entries) && return (0, nothing)
    oq = order(orthogonal_group(discriminant_group(genus_lattice(g))))
    patched = 0
    for (lab, im) in entries
        divides(oq, im)[1] || return (patched, "$lab: image order $im does not divide |O(q_L)| = $oq")
    end
    for (lab, im) in entries
        path = label_path(root, "lattice_basic_data", lab)
        isfile(path) || continue
        lines = split(read(path, String), "\n")
        fields = String.(split(lines[1], "|"))
        length(fields) == length(lfmt) || return (patched, "$path: $(length(fields)) fields, lat_basic.format has $(length(lfmt))")
        value = string(div(oq, im))
        if fields[fvi] != value
            fields[fvi] = value
            lines[1] = join(fields, "|")
            write_atomically(path, join(lines, "\n"))
            patched += 1
        end
    end
    return (patched, nothing)
end

# The genera with a sidecar and, with --legacy, every genus with a stored lattice.
function all_labels(root::String, legacy::Bool)
    labels = Set{String}()
    for (dir, strip_class) in [("festi_veniani", false), ("lattice_basic_data", true)]
        (strip_class && !legacy) && continue
        base = joinpath(root, dir)
        isdir(base) || continue
        for (_, _, files) in walkdir(base), f in files
            endswith(f, ".tmp") && continue
            push!(labels, strip_class ? f[1:findlast('.', f) - 1] : f)
        end
    end
    return sort!(collect(labels))
end

function main(args::Vector{String})
    legacy = "--legacy" in args
    rest = filter(!=("--legacy"), args)
    root = isempty(rest) ? "." : rest[1]
    gfmt = read_format(joinpath(root, "genera_basic.format"))
    lfmt = read_format(joinpath(root, "lat_basic.format"))
    fvi = findfirst(==("festi_veniani_index"), lfmt)
    labels = length(rest) > 1 ? String.(vcat([filter(!isempty, strip.(readlines(f))) for f in rest[2:end]]...)) :
                                all_labels(root, legacy)
    patched = 0
    problems = String[]
    for label in labels
        n, problem = process_genus(root, gfmt, lfmt, fvi, label; legacy)
        patched += n
        problem === nothing || push!(problems, problem)
    end
    println("festi_veniani_index: $(length(labels)) genera, $patched lattice records updated, ",
            "$(length(problems)) problems")
    foreach(p -> println("  ", p), problems)
    return isempty(problems) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
