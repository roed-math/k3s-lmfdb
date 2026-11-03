"""
LMFDBLattices: the Oscar side of the LMFDB lattice pipeline.

Holds the pipeline's Oscar computations (genus enumeration for now; see
enumerate.jl) and the Distributed worker pool that runs them (pool.jl).  It is a
package rather than a script so that its precompile workload (workload.jl) caches
the compiled code of the Oscar functions it uses: without it the first definite
`representatives` call in a fresh process spent 110 s compiling (load ~160 on
mimir), with it 0.9 s.
"""
module LMFDBLattices

using Oscar

# shared with the scripts in oscar/ (bench_oscar.jl, driver.jl, resolve_indefinite.jl)
include(joinpath(@__DIR__, "..", "..", "invariants.jl"))
include(joinpath(@__DIR__, "..", "..", "load_local.jl"))
include("enumerate.jl")
include("pool.jl")
include("workload.jl")

export enumerate_record, enumerate_task, run_pool, main_pool

end
