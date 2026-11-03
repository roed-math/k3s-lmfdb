# Test the worker pool's failure handling: julia --project=oscar oscar/test_pool.jl
using LMFDBLattices, Test
out = mktempdir()
labels = ["selftest.ok.1", "selftest.sleep.1", "selftest.sleep.60", "selftest.ok.2",
          "selftest.error.1", "selftest.alloc.6", "selftest.ok.3", "selftest.ok.4"]
counts = run_pool("selftest", labels; out_root = out, jobs = 2, timeout = 10, memlimit_gb = 4, batch = 3,
                  ledger = joinpath(out, "ledger.tsv"))
st(l) = LMFDBLattices.record_status(LMFDBLattices.out_path(out, l))
@testset "pool failure handling" begin
    @test st("selftest.sleep.60") == "timeout"      # its worker was killed after 10 s
    @test st("selftest.alloc.6") == "memory"        # over the 4 GB cap
    @test st("selftest.error.1") == "error"         # an exception is a record, not a dead worker
    @test all(st(l) == "ok" for l in labels if occursin(".ok.", l) || l == "selftest.sleep.1")
    @test counts == Dict("ok" => 5, "timeout" => 1, "memory" => 1, "error" => 1)
    # resume: nothing left to do, and nothing is rerun
    @test run_pool("selftest", labels; out_root = out, jobs = 1) == counts
end
