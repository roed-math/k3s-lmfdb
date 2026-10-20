// LabelPath centralises the data directory layout: folder/rank/nplus/label,
// where the label has the form rank.nplus.det.... (see create_genus_label).
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

AssertEqual(~results, LabelPath("lattice_basic_data", "3.3.1.1"),
    "lattice_basic_data/3/3/3.3.1.1", "lattice_basic_data, definite rank 3");
AssertEqual(~results, LabelPath("shortest", "8.8.1.1.2"),
    "shortest/8/8/8.8.1.1.2", "shortest, definite rank 8");
AssertEqual(~results, LabelPath("voronoi", "4.3.25.1"),
    "voronoi/4/3/4.3.25.1", "voronoi, indefinite rank 4 nplus 3");
AssertEqual(~results, LabelPath("genera_basic", "5.2.122.66"),
    "genera_basic/5/2/5.2.122.66", "genera_basic, indefinite rank 5 nplus 2");

Report(~results, "test_label_path");
