# Genus enumeration with Oscar, for the genera where it beats Magma: indefinite
# genera of rank >= 3 (Magma's spinor kernel overcounts some of them, see CLAUDE.md)
# and definite genera of small rank.  FillGenus reads the record this writes and
# takes the classes from it instead of enumerating; everything else per lattice is
# still computed there.
#
# Record (oscar_genus/<rank>/<nplus>/<label>):
#   line 1   status=ok classes=<h> seconds=<s> method=<definite|spinor> mass=<match|mismatch|->
#            spinor_generators=<{p,...}|\N>
#   line 1+i <gram {a,b,...} row-major>|<#Aut or \N>|<root system or \N>     (class i)
# or, when the genus could not be enumerated, a single line
#   status=<error|timeout|memory|crashed> seconds=<s> [error=<message>]
# Indefinite classes are in label order (see enumerate_indefinite); definite ones in
# Oscar's order, since FillGenus sorts them (cmp_lat).

const ENUM_FORMAT_VERSION = 1

function gram_string(L::ZZLat)
    G = gram_matrix(L)
    n = nrows(G)
    return "{" * join([string(ZZ(G[i, j])) for i in 1:n for j in 1:n], ",") * "}"     # row-major
end

# genera_basic stores the mass as "{num, den}" (or \N for indefinite genera)
function parse_stored_mass(s::AbstractString)
    s in ("\\N", "None", "") && return nothing
    a, b = parse.(BigInt, split(strip(s, ['{', '}']), ","))
    return QQ(a) // QQ(b)
end

"""
    enumerate_definite(L, stored_mass) -> (classes, mass_status)

All classes of the genus of the positive definite lattice L, each with its
automorphism group order and root system.  Oscar's `representatives` verifies the
list against the mass formula; this also compares that mass with the one the
listing stage (Sage) stored.
"""
function enumerate_definite(L::ZZLat, stored_mass)
    G = genus(L)
    reps = representatives(G)
    m = mass(G)
    tally = sum(1 // automorphism_group_order(R) for R in reps; init = QQ(0))
    tally == m || error("class list fails the mass formula: $(tally) != $(m)")
    mass_status = stored_mass === nothing ? "-" : (stored_mass == m ? "match" : "mismatch")
    classes = [(gram_string(lll(R)), string(automorphism_group_order(R)), root_string(R)) for R in reps]
    return classes, mass_status
end

"""
    enumerate_indefinite(L) -> (classes, generators)

The classes of the genus of the indefinite lattice L (rank >= 3), in label order:
by Eichler they are the improper spinor genera, built by Hecke from L with
`spinor_genera_in_genus`.  Class 1 is L itself (FillGenus passes the published rep)
and the class with spinor coordinates e in {0,1}^k, relative to the generators
q_1 < ... < q_k of `improper_spinor_generators` (the smallest primes outside the bad
primes, each kept when independent of those before it), is at position
1 + sum e_i 2^(i-1) -- the labelling rule recorded in CLAUDE.md, which Hecke's
construction order follows.  Unlike Magma's SpinorKernel, Hecke's spinor group takes
the signature into account; the count is right once Hecke.jl PR #2611 is in, which
hecke_fixes.jl applies at load time to Hecke versions without it.
"""
function enumerate_indefinite(L::ZZLat)
    gens = Hecke.improper_spinor_generators(genus(L))
    S = Hecke.spinor_genera_in_genus(L)
    length(S) == 2^length(gens) || error("$(length(S)) spinor genera for $(length(gens)) generators")
    S[1] === L || error("spinor_genera_in_genus did not start from the given lattice")
    classes = [(gram_string(R), "\\N", "\\N") for R in S]
    return classes, gens
end

"""
    enumerate_record(rank, nplus, rep, stored_mass) -> String

The oscar_genus record for the genus of the lattice with Gram matrix `rep`
(row-major entries).
"""
function enumerate_record(rank::Int, nplus::Int, rep::Vector{BigInt}, stored_mass::AbstractString)
    t0 = time()
    L = integer_lattice(; gram = matrix(ZZ, rank, rank, rep), cached = false)
    if nplus == rank
        classes, mass_status = enumerate_definite(L, parse_stored_mass(stored_mass))
        method, gens = "definite", "\\N"
    else
        rank >= 3 || error("indefinite genera of rank $rank are not enumerated with Oscar")
        classes, g = enumerate_indefinite(L)
        method, mass_status, gens = "spinor", "-", "{" * join(string.(g), ",") * "}"
    end
    head = "status=ok classes=$(length(classes)) seconds=$(round(time() - t0; digits = 3)) " *
           "method=$method mass=$mass_status spinor_generators=$gens version=$ENUM_FORMAT_VERSION"
    return join([head; [join(c, "|") for c in classes]], "\n") * "\n"
end

failure_record(status::AbstractString, seconds; message = nothing) =
    "status=$status seconds=$(round(seconds; digits = 3))" *
    (message === nothing ? "" : " error=" * replace(first(string(message), 300), r"[\s|]+" => "_")) * "\n"

"""
    enumerate_task(label, root) -> String

The pool task: read the genus's genera_basic record under `root` and enumerate it.
"""
function enumerate_task(label::AbstractString, root::AbstractString)
    fmt = read_format(joinpath(root, "genera_basic.format"))
    pieces = split(first(eachline(label_path(String(root), "genera_basic", label))), "|")
    length(pieces) == length(fmt) || error("$label: $(length(pieces)) fields, genera_basic.format has $(length(fmt))")
    rec = Dict(zip(fmt, pieces))
    rank, nplus = parse(Int, rec["rank"]), parse(Int, rec["nplus"])
    rep = parse.(BigInt, split(strip(rec["rep"], ['{', '}']), ","))
    return enumerate_record(rank, nplus, rep, rec["mass"])
end
