# One embedding job: run Stevell's overlattice / K3-embedding / corank-1 code on
# locally loaded genera.  The code lives in LMFDBLattices/src/embed.jl (with the
# pipeline's changes to the vendored code/lmfdb_lattice_enum.jl); run_all.py runs
# the jobs on the Oscar worker pool (oscar/pool.jl --task embed) instead of this
# script, which is for running a single job by hand.
#
# Usage (run from the lattices directory):
#   julia --project=oscar oscar/driver.jl --task overlat --labels-file F [options]
#   julia --project=oscar oscar/driver.jl --task k3      --labels-file F [options]
#   julia --project=oscar oscar/driver.jl --task corank1 --labels-file F --labels-file2 F2 [options]
# where F lists domain genus labels (one per line) and, for corank1, F2 lists
# the rank+1 codomain genus pool.  Options: --data-root DIR (default "."),
# --out-root DIR (default "inclusions").
#
# Resume: a job whose committed directory already holds a manifest recording the
# same inputs (sha1 of its label files) is skipped, so rerunning the embedding
# stage only computes the jobs that did not finish.  run_all.py splits each task's
# domain genera into det-sorted chunks (--embed-chunk), one job each, so a timeout
# or OOM loses one chunk rather than a whole rank or rank pair.
#
# Output layout (LabelPath sharding, keyed by the DOMAIN lattice label):
#   <out-root>/<rank>/<nplus>/<domain_label>/{prim_emb.txt, overlat.txt, k3_embed.txt}
# The line formats inside are Stevell's; write_table.jl assembles the
# lat_inclusions table from them.

using LMFDBLattices
exit(embed_cli(ARGS))
