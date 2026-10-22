# One embedding job: run Stevell's overlattice / K3-embedding / corank-1 code
# (code/lmfdb_lattice_enum.jl, vendored unmodified) on locally loaded genera.
#
# Usage (run from the lattices directory, one invocation per GNU parallel job):
#   julia --project=oscar oscar/driver.jl --task overlat --labels-file F [options]
#   julia --project=oscar oscar/driver.jl --task k3      --labels-file F [options]
#   julia --project=oscar oscar/driver.jl --task corank1 --labels-file F --labels-file2 F2 [options]
# where F lists domain genus labels (one per line) and, for corank1, F2 lists
# the rank+1 codomain genus pool.  Options: --data-root DIR (default "."),
# --out-root DIR (default "inclusions").
#
# Resume: a job whose committed directory already holds a manifest recording the
# same inputs (sha1 of its label files) is skipped, so rerunning the embedding
# stage only computes the jobs that did not finish.  run_all.py splits each task's
# domain genera into det-sorted chunks (--embed-chunk), one job each, so a timeout
# or OOM loses one chunk rather than a whole rank or rank pair.
#
# Output layout (LabelPath sharding, keyed by the DOMAIN lattice label):
#   <out-root>/<rank>/<nplus>/<domain_label>/{prim_emb.txt, overlat.txt, k3_embed.txt}
# The line formats inside are Stevell's; write_table.jl assembles the
# lat_inclusions table from them.  The save_data_lmfdb methods below override the
# vendored ones (same signatures, defined after the include) to fix two issues:
# the rank-in-filename parse that breaks for rank >= 10, and the flat data/rankN
# layout that puts every label in one directory.

using SHA
include(joinpath(@__DIR__, "load_local.jl"))
include(joinpath(@__DIR__, "code", "lmfdb_lattice_enum.jl"))

# ---- overrides of the vendored writers (include order makes these win) ----

function _inclusion_file(out_root::String, domain_label::AbstractString, fname::String)
    parts = split(domain_label, ".")
    dir = joinpath(out_root, String(parts[1]), String(parts[2]), String(domain_label))
    mkpath(dir)
    return joinpath(dir, fname)
end

# corank-1 primitive embedding: x = (codomain label, domain label, complement coords)
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, ZZMatrix},
)
  _path = _inclusion_file(_save_path, x[2], "prim_emb.txt")
  str = x[1] * "_["
  v = x[3]
  for j in 1:ncols(v)
    str *= "$(v[1,j])"
    if j != ncols(v)
      str *= ","
    end
  end
  str *= "]\n"
  open(_path, "a") do f
    write(f, str)
  end
end

# prime-index overlattice: x = (domain label, codomain label, p, glue vector)
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, Hecke.IntegerUnion, Vector{FqFieldElem}},
)
  _path = _inclusion_file(_save_path, x[1], "overlat.txt")
  str = x[2] * "_$(x[3])_["
  v = x[4]
  for j in 1:length(v)
    str *= "$(v[j])"
    if j != length(v)
      str *= ","
    end
  end
  str *= "]\n"
  open(_path, "a") do f
    write(f, str)
  end
end

# K3 embedding count: x = (domain label, count)
function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, Int},
)
  _path = _inclusion_file(_save_path, x[1], "k3_embed.txt")
  open(_path, "w") do f
    write(f, "$(x[2])")
  end
end

# Override of the vendored per-genus corank-1 dispatcher (same signature, defined
# after the include, so this method wins).  Two changes:
#
#  * Dispatch works around an Oscar bug: the unimodular fast path silently loses
#    even sublattices of ODD unimodular lattices -- primitive_embeddings(I2, <2>)
#    returns (false, []) even though Z(1,1) is such an embedding -- while the
#    generic machinery handles the same pair correctly (and the even->odd coprime
#    and generic strata check out on known-true probes).  So even G2 into odd
#    unimodular G1 goes through _primitive_embeddings_generic_safe.  The ok flag
#    itself is just "found any" (the _safe functions end with
#    `return !isempty(results), results`), so ok = false needs no logging.
#
#  * A failed isometry match raises a labelled error instead of hitting an
#    undefined variable.
function _primitive_embeddings_lmfdb(
  G1s::Vector{ZZGenus},
  G2::ZZGenus,
  _save_path::String,
)
  G1s = filter(>(rank(G2))∘rank, G1s)
  filter!(G1 -> reduce(&, signature_pair(G1).>=signature_pair(G2)), G1s)
  if !is_even(G2)
    filter!(!is_even, G1s)
  end
  Ctx = Oscar.ZZLatGluingCtx()
  Ms = get_attribute(G2, :representatives)
  glabel(G) = join(split(String(get_attribute(first(get_attribute(G, :representatives)), :lmfdb_label)), ".")[1:end-1], ".")
  for G1 in G1s
    if is_unimodular(G1) && !(is_even(G2) && !is_even(G1))
      ok, res = Oscar._primitive_embeddings_in_unimodular_safe(G1, G2, Ms)
    elseif isone(gcd(numerator(det(G1)), numerator(det(G2)))) && !is_unimodular(G1)
      ok, res = Oscar._primitive_embeddings_coprime_det_safe(G1, G2, Ms)
    else
      ok, res = Oscar._primitive_embeddings_generic_safe(G1, G2, Ms; Ctx, vi=(1, -1))
    end
    isempty(res) && continue
    Ls = get_attribute(G1, :representatives)
    for x in res
      @assert genus(x[1]) == G1
      j = findfirst(isequal(gram_matrix(x[2]))∘gram_matrix, Ms)
      @assert !isnothing(j)
      M = Ms[j]
      label_bottom = get_attribute(M, :lmfdb_label)
      matched = false
      for L in Ls
        if is_definite(G1)
          ok2, f_top = is_isometric_with_isometry(x[1], L)
        elseif length(Ls) == 1
          ok2 = true
        else
          ok2 = is_isometric(x[1], L)
        end
        !ok2 && continue
        matched = true
        label_top = get_attribute(L, :lmfdb_label)
        if !is_definite(G1)
          v = zero_matrix(ZZ, 0, 0)
        else
          v = map_entries(ZZ, coordinates(basis_matrix(x[3]), x[1])*f_top)
        end
        save_data_lmfdb(_save_path, (String(label_top), String(label_bottom), v))
      end
      matched || error("no class of ", glabel(G1), " isometric to an embedding codomain for domain ", label_bottom)
    end
  end
  return nothing
end

# Copy of the vendored _prime_index_overlattices_lmfdb that iterates over the
# domain genera Gs but looks overlattices up in `pool` (every genus of the rank).
# The vendored routine searches Gs itself, which is right only when Gs is the whole
# rank: an overlattice has the same rank and a smaller det, so with det-sorted
# chunks its genus is usually in an earlier chunk, and the row was dropped (983 of
# 32710 rows at R=4, D=120).  Only the lookups (findfirst, pool[k]) differ from the
# vendored code.
function _prime_index_overlattices_pool(
  Gs::Vector{ZZGenus},
  pool::Vector{ZZGenus},
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
            k = findfirst(isequal(G2), pool)
            isnothing(k) && continue
            # Recover the LMFDB label of L
            if length(get_attribute(pool[k], :representatives)) == 1
              label_top = get_attribute(only(get_attribute(pool[k], :representatives)), :lmfdb_label)
            else
              for L2 in get_attribute(pool[k], :representatives)
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

# ---- argument parsing ----

function parse_args(args)
    opts = Dict{String, String}("data-root" => ".", "out-root" => "inclusions")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") || error("unexpected argument $(args[i])")
        key = args[i][3:end]
        i + 1 <= length(args) || error("missing value for --$key")
        opts[key] = args[i + 1]
        i += 2
    end
    return opts
end

function inputs_digest(opts)
    files = [opts["labels-file"]]
    haskey(opts, "labels-file2") && opts["labels-file2"] != "-" && push!(files, opts["labels-file2"])
    return bytes2hex(sha1(join([read(f, String) for f in files], "\0")))
end

# The first manifest line is "key=value key=value ..."; return the value of key.
function manifest_field(path, key)
    isfile(path) || return nothing
    for kv in split(first(eachline(path)))
        k, _, v = partition_kv(kv)
        k == key && return v
    end
    return nothing
end
partition_kv(kv) = (i = findfirst('=', kv); isnothing(i) ? (kv, "", "") : (kv[1:i-1], "=", kv[i+1:end]))

function main()
    t0 = time()
    opts = parse_args(ARGS)
    task = opts["task"]
    root = opts["data-root"]
    # Atomic staging: everything is written under staging/<jobtag>/ and moved to
    # data/<jobtag>/ in one rename only after the whole job (manifest included)
    # finished.  A killed or failed job leaves only an ignored staging directory,
    # and a retry re-stages from scratch, so partial output is never consumed and
    # a rerun can never append duplicate orbits into the multiplicity files.
    jobtag = string(task, "_", replace(basename(opts["labels-file"]), ".txt" => ""))
    stage = joinpath(opts["out-root"], "staging", jobtag)
    committed = joinpath(opts["out-root"], "data", jobtag)
    digest = inputs_digest(opts)
    if manifest_field(joinpath(committed, "manifest.txt"), "inputs") == digest
        println("SKIP task=$task job=$jobtag (already committed with the same inputs)")
        return
    end
    ispath(stage) && rm(stage; recursive = true)
    mkpath(stage)
    out = stage
    labels = filter(!isempty, strip.(readlines(opts["labels-file"])))
    has_pool = get(opts, "labels-file2", "-") != "-"
    if task == "overlat" && has_pool
        # Load the rank's pool once and take this chunk's domains from it.
        pool_labels = filter(!isempty, strip.(readlines(opts["labels-file2"])))
        pool, pool_skipped = load_genera(root, String.(pool_labels))
        mine = Set(labels)
        glab(G) = join(split(String(get_attribute(first(get_attribute(G, :representatives)), :lmfdb_label)), ".")[1:end-1], ".")
        Gs = filter(G -> glab(G) in mine, pool)
        # skipped: this chunk's unloadable domains plus the pool's other ones
        # (overlattices into those genera cannot be labelled)
        skipped = pool_skipped
    else
        Gs, skipped = load_genera(root, String.(labels))
    end
    n_skipped = length(skipped)
    if task == "overlat"
        isempty(Gs) || (has_pool ? _prime_index_overlattices_pool(Gs, pool, out)
                                 : _prime_index_overlattices_lmfdb(Gs, out))
    elseif task == "k3"
        isempty(Gs) || _embeddings_in_K3_lattice(Gs, out)
    elseif task == "corank1"
        labels2 = filter(!isempty, strip.(readlines(opts["labels-file2"])))
        G1s, skipped2 = load_genera(root, String.(labels2))
        n_skipped += length(skipped2)
        if !isempty(Gs) && !isempty(G1s)
            # loop the per-genus method ourselves: the vendored Vector{ZZGenus}
            # wrapper append!s the per-genus return value, which is `nothing`
            # (data is written through save_data_lmfdb), and crashes
            for G2 in Gs
                _primitive_embeddings_lmfdb(G1s, G2, out)
                GC.gc()
            end
        end
    else
        error("unknown task $task (expected overlat, k3 or corank1)")
    end
    # Completion manifest, staged with the data: its presence in data/<jobtag>/
    # certifies that every loaded genus of this job was fully processed.  The
    # loader's skipped labels are recorded so the certificate derivation
    # (oscar/refine_universality.py) can cap the certified det range below the
    # smallest skipped genus.
    open(joinpath(stage, "manifest.txt"), "w") do f
        # (fields after skipped= are informational; readers take lines 2.. as labels)
        write(f, "task=$task loaded=$(length(Gs)) skipped=$(n_skipped) inputs=$digest seconds=$(round(time() - t0; digits = 1))\n")
        for s in skipped
            write(f, s, "\n")
        end
        if task == "corank1"
            for s in skipped2
                write(f, s, "\n")
            end
        end
    end
    mkpath(dirname(committed))
    ispath(committed) && rm(committed; recursive = true)
    mv(stage, committed)
    println("DONE task=$task genera=$(length(Gs)) skipped=$(n_skipped)")
end

main()
