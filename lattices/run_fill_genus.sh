#!/bin/bash
# Wrapper around magma run_fill_genus.m that reliably detects failures.
# Magma segfaults can report exit code 0 over SSH, so we verify that
# every label that was STARTed also got a DONE or ERROR line.
# Labels that were started but got neither are reported as segfaults.
#
# Usage:
#   ./run_fill_genus.sh LABELS [VERBOSE] [TIMEOUT]
#   where LABELS is a single label or colon-separated batch
#
# Usage with GNU Parallel (use --resume-failed, not --resume):
#   parallel --sshloginfile servers.txt --joblog jobs.log --eta --resume-failed \
#     'cd ~/projects/k3s-lmfdb/lattices && ./run_fill_genus.sh {}' \
#     :::: labels.txt

set -u

labels="$1"
verbose="${2:-0}"
timeout="${3:-60}"

cd "$(dirname "$0")" || exit 1

tmpout=$(mktemp)
trap 'rm -f "$tmpout"' EXIT

magma -b labels:="$labels" verbose:="$verbose" timeout:="$timeout" done:=1 run_fill_genus.m 2>&1 | tee "$tmpout" | grep -Ev '^(DONE|START): '
ret=${PIPESTATUS[0]}

# Find labels that were started but never completed (segfaults)
started=$(grep '^START: ' "$tmpout" | sed 's/^START: //' | sort)
finished=$(grep -E '^(DONE|ERROR): ' "$tmpout" | sed 's/^[^:]*: //' | sed 's/: .*//' | sort)
crashed=$(comm -23 <(echo "$started") <(echo "$finished"))

if [ -n "$crashed" ]; then
    while IFS= read -r l; do
        echo "ERROR: $l: segfault (no DONE/ERROR sentinel)"
    done <<< "$crashed"
    exit 1
fi

if [ $ret -ne 0 ]; then
    exit $ret
fi

exit 0
