// OrthogonalFactorLookups (connect_genus.m): which orthogonal summands connect leaves
// to the Oscar resolver (oscar/LMFDBLattices/src/resolve.jl), and how.
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();
// nothing deferred: no lookups, connect stores the factor list itself
AssertEqual(~results, OrthogonalFactorLookups(["1.1.1.3.1", "2.2.11.19.2"], ["", ""]), [Parent("") | ], "all labelled");
// one deferred: every summand becomes a lookup, a labelled one with its own label as the candidate
AssertEqual(~results, OrthogonalFactorLookups(["1.1.1.3.1", "\\N"], ["", "2.2.11.19"]),
            ["1.1.1.3;1.1.1.3.1", "2.2.11.19"], "one deferred (hash overrun)");
AssertEqual(~results, OrthogonalFactorLookups(["\\N", "\\N"], ["2.2.11.19;2.2.11.19.1,2.2.11.19.2", "1.1.1.3"]),
            ["2.2.11.19;2.2.11.19.1,2.2.11.19.2", "1.1.1.3"], "all deferred, with and without candidates");
// a final null (genus not stored): the list cannot be completed, so no lookups
AssertEqual(~results, OrthogonalFactorLookups(["\\N", "\\N"], ["", "2.2.11.19"]), [Parent("") | ], "final null summand");
Report(~results, "test_orthogonal_factor_lookups");
