# Label lookups connect could not decide, resolved with Oscar's is_isometric.
#
# connect's FindLabel identifies a lattice by its genus, the per-lattice hash, and
# an isometry test among the classes sharing the hash.  It defers the lookup (null
# label, plus a line in label_lookups/<rank>/<nplus>/<source label>) when
#   * the genus is indefinite of rank >= 3 with several classes: Magma cannot tell
#     its classes (spinor genera) apart -- every class is a candidate; or
#   * hashing the lattice or the isometry test overran its time limit in Magma,
#     whose IsIsometric failed on 9-18 of ~34 benchmark cases per rank at ranks
#     10-16 where Oscar's never did -- the classes sharing the hash are the
#     candidates (all classes if the hash itself overran).
# Each line is `column|genus label|gram|candidate labels` (the last field comma
# separated, empty for "every class").  The task finds the one candidate isometric
# to the lattice and writes its label into that column of the source lattice's
# advanced record (lattice_advanced_data), atomically.
#
# Column `orthogonal_factors` is the exception: connect writes one such line per
# orthogonal summand (a summand it could label has that label as its only
# candidate), and the factor list and multiplicities are written only once every
# summand is labelled, sorted and grouped as connect does (LatSortKey: rank, nplus
# descending, det, label).  run_all.py runs this between the two connect passes,
# so the second pass derives the decomposable fields from them.

const GENUS_CLASSES = Dict{String, Union{Nothing, Vector{ZZLat}}}()   # per worker
const FORMATS = Ref{Union{Nothing, PipelineFormats}}(nothing)

function genus_classes(root::String, genus::AbstractString)
    get!(GENUS_CLASSES, String(genus)) do
        FORMATS[] === nothing && (FORMATS[] = pipeline_formats(root))
        meta = genus_meta(root, FORMATS[], genus)
        meta === nothing ? nothing : load_class_lattices(root, FORMATS[], meta)
    end
end

function lattice_of_gram(gram::AbstractString)
    entries = parse.(BigInt, split(strip(gram, ['{', '}']), ","))
    n = isqrt(length(entries))
    return integer_lattice(; gram = matrix(ZZ, n, n, entries), cached = false)
end

"""
    resolve_task(source_label, root) -> String

Resolve the deferred lookups of one source lattice and patch its advanced record.
The record is status=ok when every lookup was resolved or had to stay null because
its genus is not completely stored (`unresolved`), and status=error when a lattice
matched no candidate or several -- an enumeration or hashing inconsistency.
"""
function resolve_task(source::AbstractString, root::AbstractString)
    root = String(root)
    t0 = time()
    adv = read_format(joinpath(root, "lat_advanced.format"))
    path = label_path(root, "lattice_advanced_data", source)
    isfile(path) || error("$source: no advanced record")
    lines = split(read(path, String), "\n")
    fields = String.(split(lines[1], "|"))               # split keeps empty fields
    length(fields) == length(adv) || error("$path: $(length(fields)) fields, lat_advanced.format has $(length(adv))")
    resolved = unresolved = 0
    problems = String[]
    factors = String[]                  # orthogonal_factors lookups, labelled
    factor_lines = 0
    for line in eachline(label_path(root, "label_lookups", source))
        isempty(line) && continue
        parts = split(line, "|")
        col, genus, gram = parts[1], parts[2], parts[3]
        cands = length(parts) >= 4 && !isempty(parts[4]) ? Set(split(parts[4], ",")) : nothing
        ci = findfirst(==(col), adv)
        ci === nothing && error("$source: unknown column $col")
        col == "orthogonal_factors" && (factor_lines += 1)
        Ls = genus_classes(root, genus)
        if Ls === nothing                                  # not all classes stored
            unresolved += 1
            continue
        end
        Q = lattice_of_gram(gram)
        pool = cands === nothing ? Ls : [L for L in Ls if get_attribute(L, :lmfdb_label) in cands]
        hits = [L for L in pool if is_isometric(Q, L)]
        if length(hits) != 1
            push!(problems, "$col:$(length(hits))_of_$(length(pool))_candidates_in_$genus")
            continue
        end
        if col == "orthogonal_factors"
            push!(factors, get_attribute(hits[1], :lmfdb_label))
            continue
        end
        fields[ci] = get_attribute(hits[1], :lmfdb_label)
        resolved += 1
    end
    if factor_lines > 0 && length(factors) == factor_lines
        fields[findfirst(==("orthogonal_factors"), adv)], fields[findfirst(==("orthogonal_multiplicities"), adv)] =
            factor_columns(factors)
        resolved += factor_lines
    end
    if resolved > 0
        lines[1] = join(fields, "|")
        write_record(path, join(lines, "\n"))
    end
    status = isempty(problems) ? "ok" : "error"
    return "status=$status seconds=$(round(time() - t0; digits = 3)) resolved=$resolved unresolved=$unresolved" *
           (isempty(problems) ? "" : " error=" * join(problems, ";")) * "\n"
end

# The orthogonal_factors and orthogonal_multiplicities values for these summand
# labels, as connect writes them: CountFibers by LatSortKey, i.e. the distinct labels
# sorted by (rank, -nplus, det, label) with their counts.
function factor_columns(labels::Vector{String})
    key(l) = (p = split(l, "."); (parse(Int, p[1]), -parse(Int, p[2]), parse(BigInt, p[3]), l))
    distinct = sort(unique(labels); by = key)
    return "{" * join(distinct, ",") * "}", "{" * join([count(==(l), labels) for l in distinct], ",") * "}"
end
