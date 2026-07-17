# Local loader replacing LMFDB.jl: builds the Vector{ZZGenus} inputs that
# code/lmfdb_lattice_enum.jl expects -- each genus carrying the attribute
# :representatives (all classes, as ZZLat) and each class the attribute
# :lmfdb_label -- straight from the pipeline's flat files, so no database
# connection (and no dev LMFDB.jl) is needed.
#
# A genus qualifies only if it is COMPLETELY enumerated: its genera_advanced
# record exists with class_number >= 1 and every lattice_basic_data file
# genus_label.1 .. genus_label.class_number is present.  Feeding partial class
# lists to the embedding code would silently produce wrong or missing rows.

using Oscar

struct GenusMeta
    label::String
    rank::Int
    nplus::Int
    class_number::Int
end

# <folder>/<rank>/<nplus>/<identifier>, the pipeline's LabelPath convention
function label_path(root::String, folder::String, label::AbstractString)
    parts = split(label, ".")
    return joinpath(root, folder, String(parts[1]), String(parts[2]), String(label))
end

read_format(path::String) = String.(split(strip(first(eachline(path))), "|"))

struct PipelineFormats
    basic::Vector{String}                  # genera_basic.format
    combined::Vector{String}               # genera_advanced files hold basic + advanced fields
    lat::Vector{String}                    # lat_basic.format
end

function pipeline_formats(root::String)
    basic = read_format(joinpath(root, "genera_basic.format"))
    adv = read_format(joinpath(root, "genera_advanced.format"))
    lat = read_format(joinpath(root, "lat_basic.format"))
    return PipelineFormats(basic, vcat(basic, adv), lat)
end

function genus_meta(root::String, fmts::PipelineFormats, label::AbstractString)
    apath = label_path(root, "genera_advanced", label)
    isfile(apath) || return nothing
    pieces = split(first(eachline(apath)), "|")
    length(pieces) == length(fmts.combined) || return nothing
    fields = Dict(zip(fmts.combined, pieces))
    cn = get(fields, "class_number", "\\N")
    cn in ("\\N", "None", "") && return nothing
    return GenusMeta(String(label), parse(Int, fields["rank"]), parse(Int, fields["nplus"]),
                     parse(Int, cn))
end

function load_class_lattices(root::String, fmts::PipelineFormats, meta::GenusMeta)
    grami = findfirst(==("gram"), fmts.lat)
    Ls = ZZLat[]
    for i in 1:meta.class_number
        lab = string(meta.label, ".", i)
        p = label_path(root, "lattice_basic_data", lab)
        isfile(p) || return nothing                        # incomplete enumeration
        pieces = split(first(eachline(p)), "|")
        entries = parse.(BigInt, split(strip(pieces[grami], ['{', '}']), ","))
        n = meta.rank
        length(entries) == n^2 || return nothing
        L = integer_lattice(gram = matrix(ZZ, n, n, entries))
        set_attribute!(L, :lmfdb_label => lab)
        push!(Ls, L)
    end
    return Ls
end

"""
    load_genera(root, labels) -> (Vector{ZZGenus}, skipped::Vector{String})

Build attribute-carrying genera for the given genus labels.  Labels whose
enumeration is missing or incomplete are returned in `skipped` (they must be
left out rather than half-loaded).
"""
function load_genera(root::String, labels::Vector{String})
    fmts = pipeline_formats(root)
    Gs = ZZGenus[]
    skipped = String[]
    for label in labels
        isempty(label) && continue
        meta = genus_meta(root, fmts, label)
        if meta === nothing || meta.class_number < 1
            push!(skipped, label)
            continue
        end
        Ls = load_class_lattices(root, fmts, meta)
        if Ls === nothing
            push!(skipped, label)
            continue
        end
        G = genus(Ls[1])
        set_attribute!(G, :representatives => Ls)
        push!(Gs, G)
    end
    return Gs, skipped
end
