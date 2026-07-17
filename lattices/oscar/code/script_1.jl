using Oscar, LMFDB # Those are DEV
db = LMFDB.LMFDBLite.LMFDBConnection();
l1 = LMFDB.genera(db; rank=1, det = <=(500));
filter!(Base.Fix2(has_attribute, :representatives), l1);
l2 = LMFDB.genera(db; rank=2, det = >=(-250) & <=(250));
filter!(Base.Fix2(has_attribute, :representatives), l2);
include("path/to/file/lmfdb_lattice_enum.jl")
_compute_overlat_and_k3_embed(l1)
_compute_overlat_and_k3_embed(l2)
_compute_prim_embed(l2, l1) 
