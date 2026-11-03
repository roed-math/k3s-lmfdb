"""
LMFDBLattices: the Oscar side of the LMFDB lattice pipeline.

Holds the pipeline's Oscar computations (genus enumeration, enumerate.jl; label
lookups connect could not decide, resolve.jl; the embedding stage, embed.jl) and the Distributed worker pool that runs them (pool.jl).
Loading it also fixes a Hecke bug in the spinor genera of indefinite genera
(hecke_fixes.jl, Hecke.jl PR #2611).  It is a
package rather than a script so that its precompile workload (workload.jl) caches
the compiled code of the Oscar functions it uses: without it the first definite
`representatives` call in a fresh process spent 110 s compiling (load ~160 on
mimir), with it 0.9 s.
"""
module LMFDBLattices

using Oscar
using SHA

# shared with the scripts in oscar/ (bench_oscar.jl, driver.jl)
include(joinpath(@__DIR__, "..", "..", "invariants.jl"))
include(joinpath(@__DIR__, "..", "..", "load_local.jl"))
include("hecke_fixes.jl")
include("enumerate.jl")
include("resolve.jl")
include("embed.jl")
include("pool.jl")
include("workload.jl")

# Hecke fixes cannot be made while this package precompiles (no method overwriting
# then), so they are made on every load, in the master and in each pool worker.
function __init__()
    apply_hecke_pr2611!()
end

export enumerate_record, enumerate_task, resolve_task, embed_task, embed_job, embed_cli, run_pool, main_pool

end
