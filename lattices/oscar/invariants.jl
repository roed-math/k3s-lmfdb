# Lattice invariants computed with Oscar, as strings in the formats the Magma
# pipeline produces for the same quantities, so the two can be compared value by
# value (benchmark/bench_oscar.jl) and so these can later fill pipeline columns.
#
# Oscar has no theta series, canonical form or Hecke (adjacency) matrices for
# integer lattices; everything else here is a thin wrapper around a documented
# Oscar/Hecke function.  Most of these need a definite lattice (as do the Magma
# intrinsics they mirror): kissing_number, minimum, successive_minima,
# root_lattice_recognition, irreducible_components and the theta series below
# check it themselves.

using Oscar

lattice_from_gram(n::Int, entries::AbstractVector) =
    integer_lattice(; gram = matrix(QQ, n, n, entries), cached = false)

# A rational as the pipeline writes it: "3" or "3/4".
qq_string(x) = (x = QQ(x); isone(denominator(x)) ? string(numerator(x)) : string(numerator(x), "/", denominator(x)))

"""
    theta_series(L, prec) -> Vector{ZZRingElem}

The coefficients a_0, ..., a_prec of the theta series sum_v q^(v,v) of the positive
definite lattice L (the exponent is the norm, as in Magma's ThetaSeries(L, prec)).
short_vectors_iterator returns one vector of each pair +-v.
"""
function theta_series(L::ZZLat, prec::Int)
    @req is_positive_definite(L) "theta series need a positive definite lattice"
    c = zeros(ZZRingElem, prec + 1)
    c[1] = 1
    prec == 0 && return c
    for (_, nv) in short_vectors_iterator(L, prec, Int)
        c[Int(ZZ(nv)) + 1] += 2
    end
    return c
end

const ROOT_TYPES = (:I, :A, :D, :E)

"""
    root_string(L) -> String

The root system of the sublattice spanned by the vectors of norm 1 and 2 of the
definite lattice L, in the pipeline's root_sublattice format (RootString in
orth_dec.mag): components ordered I, A, D, E and by rank within a type, multiplicity
prefixes, all norm-1 components merged into one I_k, and "0" when there are none.
"""
function root_string(L::ZZLat)
    types, _ = root_lattice_recognition(L)
    k = 0
    comps = Tuple{Int, Int}[]
    for (t, n) in types
        if t === :I
            k += n
            continue
        end
        t === :D && n == 3 && (t = :A)              # D3 = A3; the pipeline says A3
        push!(comps, (findfirst(==(t), ROOT_TYPES), n))
    end
    k > 0 && push!(comps, (1, k))
    isempty(comps) && return "0"
    sort!(comps)
    parts = String[]
    i = 1
    while i <= length(comps)
        j = i
        while j < length(comps) && comps[j + 1] == comps[i]
            j += 1
        end
        mult = j - i + 1
        push!(parts, string(mult > 1 ? mult : "", ROOT_TYPES[comps[i][1]], comps[i][2]))
        i = j + 1
    end
    return join(parts, " ")
end

# An orthogonal decomposition as its sorted multiset of (rank, determinant) pairs;
# the Magma side summarises OrthogonalDecomposition the same way.
decomposition_string(comps) =
    join(["$(r):$(d)" for (r, d) in sort([(rank(c), ZZ(abs(det(c)))) for c in comps])], ",")

successive_minima_string(L::ZZLat) = join(qq_string.(successive_minima(L)), ",")

"""
    aut_label(G) -> Union{String, Nothing}

The SmallGroups label "order.id" of the finite group G, or nothing when GAP's small
group identification does not cover its order (the pipeline's aut_label is then null).
"""
function aut_label(G)
    o = order(G)
    has_small_group_identification(o) || return nothing
    n, k = small_group_identification(G)
    return "$(n).$(k)"
end

"""
    discriminant_image(L) -> (image order, |O(q_L)|)

The order of the image of O(L) in the orthogonal group O(q_L) of the discriminant
form, and the order of O(q_L) itself; their quotient is the Festi-Veniani index
[O(q_L) : im O(L)].  Works for indefinite L too (via strong approximation).
"""
function discriminant_image(L::ZZLat)
    imOL, inj = image_in_Oq(L)
    return order(imOL), order(codomain(inj))
end
