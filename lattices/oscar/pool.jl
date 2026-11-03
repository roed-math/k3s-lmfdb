# Run an LMFDBLattices task over a list of genus labels on a pool of persistent
# Oscar workers (see LMFDBLattices/src/pool.jl).
#
#   julia --project=oscar oscar/pool.jl --task enumerate --labels-file F \
#       [--jobs 4] [--timeout 600] [--memlimit 16] [--batch 20] \
#       [--data-root .] [--out-root oscar_genus] [--ledger FILE] [--retry-failed]
#
# Load the package (and with it Oscar) here, before the workers are started:
# addprocs installs its own SIGCHLD handling, after which a dev checkout of Oscar
# deadlocks in __init__ waiting for `git` (seen with Distributed elsewhere).
using LMFDBLattices
exit(main_pool(ARGS))
