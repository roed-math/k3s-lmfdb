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
    end
end
