#!/bin/bash
# Wrapper around magma run_fill_genus.m that reliably detects failures.
# Magma segfaults can report exit code 0 over SSH, so we verify that
# every label in the batch got a DONE or ERROR line, and report every
# label that did not:
#   - STARTed but never finished: Magma died while filling that genus;
#   - never STARTed: Magma died earlier in the batch (e.g. in the orth-cache
#     pre-warm, which runs before the per-label loop), so the genus was never
#     attempted.  Checking only STARTed labels would miss these entirely.
#
# Usage:
#   ./run_fill_genus.sh labels:=A:B:C [key:=value ...]
#       Every argument is passed to run_fill_genus.m as is (timeout:=, masslimit:=,
#       sizelimit:=, timelimit:=, ...); done:=1 is added.  This is how run_all.py
#       calls it, so the wrapper runs with exactly the guards run_all.py sets.
#   ./run_fill_genus.sh LABELS [VERBOSE] [TIMEOUT]
#       Old positional form (no fill guards beyond the timeout, default 60s).
#
# Usage with GNU Parallel (use --resume-failed, not --resume):
#   parallel --sshloginfile servers.txt --joblog jobs.log --eta --resume-failed \
#     'cd ~/projects/k3s-lmfdb/lattices && ./run_fill_genus.sh labels:={}' \
#     :::: labels.txt

set -u

if [[ "${1:-}" == *":="* ]]; then
    magma_args=("$@")
else
    magma_args=("labels:=$1" "verbose:=${2:-0}" "timeout:=${3:-60}")
fi

labels=""
for a in "${magma_args[@]}"; do
    case "$a" in
        labels:=*) labels="${a#labels:=}" ;;
        label:=*)  labels="${a#label:=}" ;;
    esac
done
if [ -z "$labels" ]; then
    echo "run_fill_genus.sh: no labels:= argument" >&2
    exit 2
fi

cd "$(dirname "$0")" || exit 1

tmpout=$(mktemp)
trap 'rm -f "$tmpout"' EXIT

magma -b "${magma_args[@]}" done:=1 run_fill_genus.m 2>&1 | tee "$tmpout" | grep -Ev '^(DONE|START): '
ret=${PIPESTATUS[0]}

expected=$(tr ':' '\n' <<< "$labels" | grep -v '^$' | sort -u)
started=$(grep '^START: ' "$tmpout" | sed 's/^START: //' | sort -u)
finished=$(grep -E '^(DONE|ERROR): ' "$tmpout" | sed 's/^[^:]*: //' | sed 's/: .*//' | sort -u)
unfinished=$(comm -23 <(echo "$expected") <(echo "$finished"))

if [ -n "$unfinished" ]; then
    while IFS= read -r l; do
        if grep -qxF -- "$l" <<< "$started"; then
            echo "ERROR: $l: magma died while filling this genus (exit $ret; no DONE/ERROR sentinel)"
        else
            echo "ERROR: $l: not reached (magma exited $ret before starting this genus)"
        fi
    done <<< "$unfinished"
    exit 1
fi

exit $ret
