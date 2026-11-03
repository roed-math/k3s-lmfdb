# Precompile workload: run every task on small fixed lattices once, so the compiled
# code lands in the package image instead of being compiled again in every worker.
using PrecompileTools

@setup_workload begin
    definite = BigInt[2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 6]      # a small rank-4 definite lattice
    scaled = BigInt[4, 2, 2, 4]                                            # 2 * A2
    two_classes = BigInt[1, 0, 0, 0, 20, 0, 0, 0, -25]                    # two spinor genera
    one_class = BigInt[1, 0, 0, 0, 1, 0, 0, 0, -3]
    @compile_workload begin
        enumerate_record(4, 4, definite, "\\N")
        enumerate_record(2, 2, scaled, "\\N")
        enumerate_record(3, 2, two_classes, "\\N")
        enumerate_record(3, 2, one_class, "\\N")
        L = integer_lattice(; gram = matrix(ZZ, 4, 4, definite), cached = false)
        irreducible_components(L)
        is_isometric(L, L)
        failure_record("error", 0.0; message = "x")
        # embedding stage: the corank-1 and overlattice machinery is what each
        # embedding job used to compile afresh (~270 s per job under load)
        try
            lab(G, l) = (for (i, L) in enumerate(get_attribute(G, :representatives)); set_attribute!(L, :lmfdb_label => "$l.$i"); end; G)
            mk(g, l) = (L = integer_lattice(; gram = matrix(ZZ, isqrt(length(g)), isqrt(length(g)), g)); G = genus(L);
                        set_attribute!(G, :representatives => [L]); lab(G, l))
            out = mktempdir()
            G4, G1 = mk([4], "1.1.4.0"), mk([1], "1.1.1.0")
            _prime_index_overlattices_pool([G4], [G4, G1], out)
            _primitive_embeddings_lmfdb([mk([1, 0, 0, 1], "2.2.1.0"), mk([1, 0, 0, 2], "2.2.2.0")], G1, out)
            rm(out; recursive = true)
        catch
        end
    end
end
