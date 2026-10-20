#!/usr/bin/env bash
#
# tests/run.sh — run the Magma test suites. Run from lattices/ (where
# lattices.spec lives), same as every other script in this directory.
#
#   tests/run.sh                 run the fast suites (the default)
#   tests/run.sh --slow          run ALL suites, fast and slow
#   tests/run.sh <name>          run one suite, streaming its output
#   tests/run.sh --help          this message
#
# --slow means "also include the slow suites", not "run only the slow ones".
#
# Each suite runs in its own Magma process, so global state cannot leak
# between suites and one crash cannot take out the run. Full transcripts go
# to tests/logs/<suite>.log.
#
# A suite's status comes from the completion marker its Report call prints,
# combined with the process exit code:
#
#   exactly 1 marker, PASS, exit 0        -> PASS
#   exactly 1 marker, FAIL, exit 1        -> FAIL
#   anything else (including 0 markers)   -> INCOMPLETE
#
# INCOMPLETE means test execution itself failed: a crash, a broken fixture, a
# suite that asserted nothing, or output redirection that swallowed the
# marker. It is deliberately distinct from FAIL, which means the code under
# test produced a wrong answer.
#
# Adapted from assaferan/AtkinLehnerQuotients's tests/run.sh (itself from
# assaferan/ShimuraCurveALQuotients); only the root-directory check and the
# SLOW_SUITES list differ.

set -uo pipefail

# CanonicalForm on the matrices in these suites is the cost. Measured
# 2026-09-03: test_brandhorst 238s, test_l12a_plus 8602s (~2h23m, of which
# 2930s is the `time cans := ...` line and the rest is the test_canonical
# loop -- note the in-file comments predating this claim 376s and 103s, which
# are the CanonicalForm timings alone, not the whole suite). The next
# slowest suite is test_genus_reps_orth at ~17s, so the gap is enormous;
# these two are why --slow exists.
SLOW_SUITES="test_brandhorst test_l12a_plus"

MAGMA="${MAGMA:-magma}"
LOGDIR="tests/logs"
MARKER='^@@TEST_RESULT@@'

usage() {
    # BSD sed has no \s, so the prefix strip is spelled out longhand.
    sed -n '3,10p' "$0" | sed 's|^# \{0,1\}||'
}

if [[ ! -f lattices.spec ]]; then
    echo "ERROR: run.sh must be run from the lattices/ directory." >&2
    echo "       (expected lattices.spec in \$(pwd)=$(pwd))" >&2
    exit 2
fi

# All suites, in glob order. tests/assertions.m and any non-test_* helper
# file (e.g. tests/_l12a_plus_data.m) do not match and are skipped.
all_suites() {
    local f found=0
    for f in tests/test_*.m; do
        [[ -e "$f" ]] || continue
        found=1
        basename "$f" .m
    done
    [[ $found -eq 1 ]]
}

is_slow() {
    local s
    for s in $SLOW_SUITES; do
        [[ "$1" == "$s" ]] && return 0
    done
    return 1
}

# Run one suite. Writes the transcript to $2, and sets RUN_STATUS and RUN_SECS.
# $3 = "stream" to also send output to the terminal as it runs.
#
# The results come back in globals rather than on stdout because stream mode
# pipes Magma through tee: if this function were called in a subshell or
# process substitution to capture an echoed result, that tee output would be
# captured too and would never reach the terminal.
RUN_STATUS=""
RUN_SECS=0

run_suite() {
    local name="$1" log="$2" mode="${3:-quiet}"
    local file="tests/$name.m"
    local start end rc markers word

    start=$(date +%s)
    if [[ "$mode" == "stream" ]]; then
        "$MAGMA" -b "$file" </dev/null 2>&1 | tee "$log"
        rc=${PIPESTATUS[0]}
    else
        "$MAGMA" -b "$file" </dev/null >"$log" 2>&1
        rc=$?
    fi
    end=$(date +%s)

    markers=$(grep -c "$MARKER" "$log")
    if [[ "$markers" -ne 1 ]]; then
        RUN_STATUS=INCOMPLETE
    else
        word=$(grep "$MARKER" "$log" | awk '{print $2}')
        if [[ "$word" == "PASS" && $rc -eq 0 ]]; then
            RUN_STATUS=PASS
        elif [[ "$word" == "FAIL" && $rc -eq 1 ]]; then
            RUN_STATUS=FAIL
        else
            RUN_STATUS=INCOMPLETE
        fi
    fi

    RUN_SECS=$((end - start))
}

# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------
include_slow=0
single=""

case "${1:-}" in
    --help|-h) usage; exit 0 ;;
    --slow)    include_slow=1 ;;
    "")        ;;
    -*)        echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    *)         single="$1" ;;
esac

mkdir -p "$LOGDIR"

# ---------------------------------------------------------------------------
# Single-suite mode: stream, do not summarise.
# ---------------------------------------------------------------------------
if [[ -n "$single" ]]; then
    single="${single%.m}"
    single="${single#tests/}"
    if [[ ! -f "tests/$single.m" ]]; then
        echo "ERROR: unknown test suite: $single" >&2
        echo "       (no such file: tests/$single.m)" >&2
        echo "Available suites:" >&2
        all_suites | sed 's/^/  /' >&2
        exit 2
    fi
    run_suite "$single" "$LOGDIR/$single.log" stream
    echo
    printf '%-11s %s   %ss\n' "$RUN_STATUS" "$single" "$RUN_SECS"
    [[ "$RUN_STATUS" == "PASS" ]] && exit 0
    exit 1
fi

# ---------------------------------------------------------------------------
# Multi-suite mode.
# ---------------------------------------------------------------------------
suites=()
while read -r s; do
    if [[ $include_slow -eq 1 ]] || ! is_slow "$s"; then
        suites+=("$s")
    fi
done < <(all_suites)

if [[ ${#suites[@]} -eq 0 ]]; then
    echo "ERROR: no test suites selected (tests/test_*.m matched nothing)." >&2
    exit 2
fi

if [[ $include_slow -eq 1 ]]; then
    echo "Running all ${#suites[@]} suites (including the slow ones)."
else
    echo "Running ${#suites[@]} fast suites. Use --slow to include the slow ones."
fi
echo

nfail=0
nincomplete=0
for s in "${suites[@]}"; do
    run_suite "$s" "$LOGDIR/$s.log"
    if [[ "$RUN_STATUS" == "PASS" ]]; then
        printf '%-11s %-34s %5ss\n' "$RUN_STATUS" "$s" "$RUN_SECS"
    else
        printf '%-11s %-34s %5ss   -> %s\n' \
            "$RUN_STATUS" "$s" "$RUN_SECS" "$LOGDIR/$s.log"
        [[ "$RUN_STATUS" == "FAIL" ]] && nfail=$((nfail + 1))
        [[ "$RUN_STATUS" == "INCOMPLETE" ]] && nincomplete=$((nincomplete + 1))
    fi
done

echo
bad=$((nfail + nincomplete))
if [[ $bad -eq 0 ]]; then
    echo "All ${#suites[@]} suites PASS."
    exit 0
fi
echo "$bad of ${#suites[@]} suites not PASS ($nfail FAIL, $nincomplete INCOMPLETE)."
exit 1
