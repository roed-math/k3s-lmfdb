# Test the Hecke fixes LMFDBLattices applies at load time (hecke_fixes.jl):
#   julia --project=oscar oscar/test_hecke_fixes.jl
using LMFDBLattices, Oscar, Test

@testset "Hecke.jl PR #2611 applied at load time" begin
    @test LMFDBLattices.HECKE_PR2611[] in (:applied, :already_in_hecke)
    # The PR's own test lattices (all from the census): two classes each, one before
    # the fix.
    for m in [[1, 0, 0, 0, -12, -102, 0, -102, -3], [-2, 4, -2, 4, 67, -17, -2, -17, 73],
              [-12, 6, 0, 6, 22, -6, 0, -6, 36], [1, 6, -3, 6, -12, 6, -3, 6, 222],
              [13, 5, 7, 5, 13, -25, 7, -25, -2], [-1, 7, 0, 7, -1, 0, 0, 0, 27],
              [1, 0, 0, 0, -12, 6, 0, 6, 213], [-2, 1, 17, 1, 13, 5, 17, 5, 13],
              [1, 15, 6, 15, 6, -9, 6, -9, 9], [-4, 14, 0, 14, -1, 0, 0, 0, 27],
              [-1, 7, 0, 7, -1, 0, 0, 0, 108], [8, 0, -14, 0, 18, -9, -14, -9, -25],
              [-1, -8, 13, -8, 11, 11, 13, 11, 50]]
        L = integer_lattice(; gram = matrix(QQ, 3, 3, m))
        @test length(genus_representatives(L)) == 2
    end
    # The pipeline's record for 3.2.10800.0001.11.01.293b: two classes, spinor generator
    # 7, as the patched Magma has it (11, which unpatched Magma picks, is automorphous:
    # the lattice is isometric to its 11-neighbours).
    rec = enumerate_record(3, 2, BigInt[1, 6, -3, 6, -12, 6, -3, 6, 222], "\\N")
    @test occursin("classes=2", rec) && occursin("spinor_generators={7}", rec)
    # a genus with one class keeps one
    @test occursin("classes=1", enumerate_record(3, 2, BigInt[1, 0, 0, 0, 1, 0, 0, 0, -3], "\\N"))
end
