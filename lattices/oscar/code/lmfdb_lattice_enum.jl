# **Need OSCAR 1.9.0-DEV
# - use dev Oscar
#
# **Need LMFDB (by thofma)
# There is a compatibility issue with the LMFDB connection, so one has to use
# the version dev version in the following way:
# - dev https://github.com/thofma/LMFDB.jl
# - in line 10 of `src/NumberField.jl`, add the signature
#   `db::LMFDB.LMFDBLite.LMFDBConnection` to the first input of the function
#   `Hecke.number_field(db, label::String).
#
# **Run the code (after the previous steps)
# - open julia, and run `using Oscar, LMFDB`
# - include the file `include("path/to/file/lmfdb_lattice_enum.jl")`
# - fetch the database: `db = LMFDB.LMFDBLite.LMFDBConnection()`
# - load some genera: z.B. `l6=LMFDB.genera(db; rank = 6, nplus=6, det = >=(1) & <=(div(300, 8, RoundUp)));`
# - call the functions needed below
#
# WARNING: the functions assume that the genera are enumerate. One
# should be careful that if it is used on a non-complete version of the
# LMFDB, some genera might not be enumerated... Those have to be filtered
# out before throwing them the algorithms (we need the LMFDB label, for
# instance).
#
# Make sure that `__save_path` is safe. See the file `writing_file.jl` to
# translate the resulting computations into a proper table to be uploaded
# on the LMFDB.

###############################################################################
#
# Main functions
#
###############################################################################

# Run the overlattice computations and the embedding into
# the K3 lattice
#
# This can be run independently for every list of genera of
# a same rank (and bounded absolute determinant)
function _compute_overlat_and_k3_embed(
  Gs::Vector{ZZGenus},
  __save_path::String = joinpath(@__DIR__, "../data"),
)
  @show __save_path
  n = rank(first(Gs))
  # Make sure that where we want to save data exists
  if !isnothing(__save_path)
    _save_path = joinpath(__save_path, "rank$n")
    isdir(_save_path) || mkdir(_save_path)
    for G in Gs, L in get_attribute(G, :representatives)
      label = get_attribute(L, :lmfdb_label)
      path_label = joinpath(_save_path, label)
      isdir(path_label) || mkdir(path_label)
    end
  else
    _save_path = nothing
  end
  _prime_index_overlattices_lmfdb(Gs, _save_path)
  _embeddings_in_K3_lattice(Gs, _save_path)
  return nothing
end

# - `G1s` consists of all genera of lattices of rank n+1, with given
#   signatures and determinant in absolute value bounded by d/(n+1)
#   for some d.
# - `G2s` consists of all genera of lattice of rank n, with same
#   possible signatures and determinant in absolute value bounded by d/n.
# - `__save_path` is the absolute path to a folder where to save data (then
#   the dispatch is made automatically; see the code below).
function _compute_prim_embed(
  G1s::Vector{ZZGenus},
  G2s::Vector{ZZGenus},
  __save_path::String = joinpath(@__DIR__, "../data"),
)
  n = rank(first(G2s))
  if !isnothing(__save_path)
    _save_path = joinpath(__save_path, "rank$n")
    isdir(_save_path) || mkdir(_save_path)
    for G2 in G2s, L in get_attribute(G2, :representatives)
      label = get_attribute(L, :lmfdb_label)
      path_label = joinpath(_save_path, label)
      isdir(path_label) || mkdir(path_label)
    end
  else
    _save_path = nothing
  end

  _primitive_embeddings_lmfdb(G1s, G2s, _save_path)
  return nothing
end

###############################################################################
#
# Saving functions
#
###############################################################################

# Save the primitive embedding data
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, ZZMatrix},
)
  m = parse(Int, first(x[1]))
  _path = joinpath(_save_path, x[2], "prim_emb_$(m).txt")
  isfile(_path) || touch(_path)
  str = x[1]*"_["
  v = x[3]
  for j in 1:ncols(v)
    k = v[1,j]
    str *= "$k"
    if j != ncols(v)
      str *= ","
    end
  end
  str *= "]\n"
  f = open(_path, "a")
  write(f, str)
  close(f)
end

# Save the overlattice data
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, Hecke.IntegerUnion, Vector{FqFieldElem}},
)
  _path = joinpath(_save_path, x[1], "overlat.txt")
  isfile(_path) || touch(_path)
  str = x[2]*"_$(x[3])_["
  v = x[4]
  for j in 1:length(v)
    k = v[j]
    str *= "$k"
    if j != length(v)
      str *= ","
    end
  end
  str *= "]\n"
  f = open(_path, "a")
  write(f, str)
  close(f)
end

# Save the K3 lattice embedding data
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, Int},
)
  _path = joinpath(_save_path, x[1], "k3_embed.txt")
  isfile(_path) || touch(_path)
  f = open(_path, "w")
  write(f, "$(x[2])")
  close(f)
end

###############################################################################
#
# Primitive embeddings
#
###############################################################################

function _primitive_embeddings_lmfdb(
  G1s::Vector{ZZGenus},
  G2s::Vector{ZZGenus},
  _save_path::String,
)
  results = Tuple{String, String, ZZMatrix}[]
  for G2 in G2s
    append!(results, _primitive_embeddings_lmfdb(G1s, G2, _save_path))
    GC.gc()
  end
  return results
end

function _primitive_embeddings_lmfdb(
  G1s::Vector{ZZGenus},
  G2::ZZGenus,
  _save_path::String,
)
  # Filter G1s to only keep genera in which G2 can embed
  # primitively
  G1s = filter(>(rank(G2))∘rank, G1s)
  filter!(G1 -> reduce(&, signature_pair(G1).>=signature_pair(G2)), G1s)
  if !is_even(G2)
    filter!(!is_even, G1s)
  end
  # Context object to avoid redundant computations
  Ctx = Oscar.ZZLatGluingCtx()

  Ms = get_attribute(G2, :representatives)
  for G1 in G1s
    # We dispatch to dedicate methods depending on the determinants of
    # G1 and G2
    if is_unimodular(G1)
      ok, res = Oscar._primitive_embeddings_in_unimodular_safe(G1, G2, Ms)
    elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
      ok, res = Oscar._primitive_embeddings_coprime_det_safe(G1, G2, Ms)
    else
      ok, res = Oscar._primitive_embeddings_generic_safe(G1, G2, Ms; Ctx, vi=(1, -1))
    end

    isempty(res) && continue
    Ls = get_attribute(G1, :representatives)
    for x in res
      @assert genus(x[1]) == G1
      # We have a primitive embedding of G2 into G1, so we need first
      # to recover which lattice of G2 we have embedded into which lattice
      # of G1
      j = findfirst(isequal(gram_matrix(x[2]))∘gram_matrix, Ms)
      @assert !isnothing(j)
      M = Ms[j]
      # LMFDB label of the domain of the embedding
      label_bottom = get_attribute(M, :lmfdb_label)
      for L in Ls
        if is_definite(G1)
          ok, f_top = is_isometric_with_isometry(x[1], L)
        elseif length(Ls) == 1
          ok = true
        else
          ok = is_isometric(x[1], L)
        end
        !ok && continue
        # LMFDB label of the codomain of the embedding
        label_top = get_attribute(L, :lmfdb_label)
        # Try to see if we can recover the coordinates of the
        # orthogonal complement of M in L after embedding
        if !is_definite(G1)
          v = zero_matrix(ZZ, 0, 0)
        else
          v = map_entries(ZZ, coordinates(basis_matrix(x[3]), x[1])*f_top)
        end
        __x = (label_top, label_bottom, v)
        save_data_lmfdb(_save_path, __x)
      end
    end
  end
  return nothing
end

###############################################################################
#
# Overlattices
#
###############################################################################

# Compute all orbits of prime index overlattices for lattices in Gs,
# which are also contained in Gs
function _prime_index_overlattices_lmfdb(
  Gs::Vector{ZZGenus},
  _save_path::String
)
  for G in Gs
    q = discriminant_group(G)
    o = order(q)
    # In those cases, the discriminant group has no
    # isotropic vectors so we can skip this part
    if isone(o)
      continue
    elseif radical(o) == o
      continue
    end
    # We want any integral overlattice, so in the
    # even case, we see the discriminant group as
    # a finite bilinear form
    if is_even(G)
      q = Hecke._as_finite_bilinear_module(q)
    end
    # Any prime index overlattice has index in pds
    pds = prime_divisors(last(elementary_divisors(q)))
    D = Dict{eltype(pds), Tuple{TorQuadModule, Vector{Tuple{TorQuadModuleMap, GAPGroupHomomorphism}}}}()
    for p in pds
      qp, _ = torsion_subgroup(q, p)
      # Function to compute 1-dimensional isotropic subspace
      # of the p-torsion part of the discriminant group
      r = Oscar._stabilizer_isotropic_elementary(qp, 1)
      D[p] = (qp, r)
    end
    # We have done the computations "locally", we can now
    # transport them to every lattice in the genus and
    # split orbits under the action of the respective
    # orthogonal groups
    for M in get_attribute(G, :representatives)
      GM, _ = image_in_Oq(M)
      # Keep track of the label of M for storage
      label_bottom = get_attribute(M, :lmfdb_label)
      T = domain(GM)
      if is_even(G)
        T = Hecke._as_finite_bilinear_module(T)
        GM = Oscar._orthogonal_group(T, matrix.(gens(GM)); check=false)
      end
      for p in pds
        K = GF(p)
        _v = zeros(K, degree(M))
        qp, tmp = D[p]
        is_empty(tmp) && continue
        Tp, TpinT = torsion_subgroup(T, p)
        OTp = orthogonal_group_bilinear(Tp)
        iso = isomorphism(PermGroup, OTp)
        Gp, _ = restrict_automorphism_group(GM, TpinT; check=false)
        isoGp, _ = iso(Gp)
        ok, phi = is_isometric_with_isometry(Tp, qp)
        @assert ok # Should always be true
        iphi = inv(phi)
        for (Hinqp, j) in tmp
          # Recover a representative for an O(q)-orbit of 1-dimensional
          # isotropic F_p-vector space
          _H, _HinTp = sub(Tp, elem_type(Tp)[phi\(qp(lift(a))) for a in gens(domain(Hinqp))])
          _, jTp = sub(codomain(iso), elem_type(codomain(iso))[iso(OTp(phi * hom(j(g)) * iphi; check=false)) for g in gens(domain(j))])
          # Split the O(q)-orbit into O(M)-orbits
          for _h in double_cosets(codomain(iso), domain(jTp), isoGp)
            h = iso\(representative(_h))
            _Hh, _ = sub(Tp, [h(Tp(lift(a))) for a in gens(_H)])
            L = cover(_Hh)
            # L is the overlattice obtained by adding a generator of
            # _Hh to M. We choose a reduced generator by multiplying
            # coordinates by p (so they are integers) and reduce them
            # modulo p (works because the basis matrix of M is the
            # identity so any vector with integral coordinates lies
            # in M)
            map!(K ∘numerator, _v, p*lift(first(gens(_Hh))))
            G2 = genus(L)
            # Try to recover L in the list G2: technically, since we
            # aim to choose Gs to contain all lattices up to a certain
            # determinant bound (in absolute value), the genus of L
            # should be in Gs. Only exception: the genus of L is not
            # yet known to be enumerated so we have already discarded it
            k = findfirst(isequal(G2), Gs)
            isnothing(k) && continue
            # Recover the LMFDB label of L
            if length(get_attribute(Gs[k], :representatives)) == 1
              label_top = get_attribute(only(get_attribute(Gs[k], :representatives)), :lmfdb_label)
            else
              for L2 in get_attribute(Gs[k], :representatives)
                !is_isometric(L, L2) && continue
                label_top = get_attribute(L2, :lmfdb_label)
              end
            end
            __x = (label_bottom, label_top, p, _v)
            save_data_lmfdb(_save_path, __x)
          end
        end
      end
    end
  end
  return nothing
end

###############################################################################
#
# Primitive embeddings in K3 lattice
#
###############################################################################

function _embeddings_in_K3_lattice(
  Gs::Vector{ZZGenus},
  _save_path::String,
)
  L = rescale(k3_lattice(), -1)
  GK3 = genus(L)
  for G in Gs
    # We can only take even lattices and we
    # are interested in lattices which are
    # not definite
    !is_even(G) && continue
    p, n = signature_pair(G)
    n in [1,2] || continue
    Ms = get_attribute(G, :representatives)
    _, r = Oscar._primitive_embeddings_in_unimodular_safe(GK3, G, Ms)
    for M in Ms
      l = count(x -> gram_matrix(x[2]) == gram_matrix(M), r)
      label_bottom = get_attribute(M, :lmfdb_label)
      __x = (label_bottom, l)
      save_data_lmfdb(_save_path, __x)
    end
  end
  return nothing
end

