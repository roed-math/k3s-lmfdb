// tests/assertions.m
//
// Shared assertion procedures. Every test suite loads this file first:
//
//     load "tests/assertions.m";
//     results := NewResults();
//     ...
//     Report(~results, "test_my_suite");
//
// Adapted from assaferan/AtkinLehnerQuotients's tests/assertions.m (itself
// from assaferan/ShimuraCurveALQuotients), with one addition: AssertTrue, for
// suites whose checks are plain boolean conditions rather than actual-vs-
// expected comparisons.
//
// SetQuitOnError(true) is set here rather than in each suite. Without it a
// Magma runtime error exits 0, so a crashed suite would look like a clean
// run; setting it here means no suite can omit it by accident. It also
// covers bare `assert` statements inside helper functions a suite calls
// (e.g. canonical_form.m's test_canonical) that predate this convention and
// were not rewritten to use AssertTrue/AssertEqual: a failure there still
// exits non-zero and prints no completion marker, so run.sh reports
// INCOMPLETE rather than a false PASS.
//
// Report emits a machine-readable completion marker immediately before
// quitting. tests/run.sh derives a suite's status from that marker combined
// with the process exit code; a missing, duplicated, or inconsistent marker
// is INCOMPLETE. Nothing else this file prints is part of that protocol, so
// the human-readable wording below can be changed freely.

SetQuitOnError(true);

TestResults := recformat< nassert : RngIntElt, failures : SeqEnum >;

function NewResults()
    return rec< TestResults | nassert := 0, failures := [] >;
end function;

// Compare one actual value against its expected value. Records a mismatch
// rather than raising, so a suite reports every failing example in one run
// instead of stopping at the first. `label` must identify the example well
// enough to locate it: include the level, discriminant, or row.
procedure AssertEqual(~results, actual, expected, label)
    results`nassert +:= 1;
    if actual ne expected then
        Append(~results`failures, <label, Sprint(expected), Sprint(actual)>);
    end if;
end procedure;

// Same bookkeeping as AssertEqual, for a plain boolean condition rather than
// an actual-vs-expected pair (e.g. "#bad_idxs eq 4", "IsIsomorphic(A, B)").
procedure AssertTrue(~results, cond, label)
    results`nassert +:= 1;
    if not cond then
        Append(~results`failures, <label, "true", "false">);
    end if;
end procedure;

// Print diagnostics, emit the completion marker, and quit with the matching
// status. Call this once, at the very end of a suite.
procedure Report(~results, name)
    // A suite that asserted nothing is a broken fixture, not a passing run.
    // Raising here emits no marker, so run.sh classifies it INCOMPLETE.
    if results`nassert eq 0 then
        error Sprintf("%o: no assertions ran", name);
    end if;

    if #results`failures eq 0 then
        printf "%o: %o assertions, all passed\n", name, results`nassert;
        printf "@@TEST_RESULT@@ PASS\n";
        quit 0;
    else
        for f in results`failures do
            printf "%o\n", f[1];
            printf "  expected: %o\n", f[2];
            printf "  actual:   %o\n", f[3];
        end for;
        printf "%o: %o assertions, %o FAILED\n",
            name, results`nassert, #results`failures;
        printf "@@TEST_RESULT@@ FAIL\n";
        quit 1;
    end if;
end procedure;
