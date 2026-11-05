// Usage: magma -b sig:=4.3 Dbound:=1000 run_tensor.m
// Here 4 is the rank of the lattices, 3 is n_+, and 1000 is the upper limit on the discriminant.  This script will raise an error if n=1 or n_+ < n/2.
// It should be run after fill_genus.m, so that lattice data is saved in 
// It will write two files to disk: tensor_indecomp/n.n+ (with the labels of all primitive indecomposable lattices) and tensor_decomp/n.n+ (with the decompositions of all primitive lattices of this signature)

AttachSpec("lattices.spec");
SetColumns(0);

n, nplus := Explode(Split(sig, "."));
n := StringToInteger(n);
nplus := StringToInteger(nplus);
Dbound := StringToInteger(Dbound);

decomp, indecomp := TensorDecompositions(<n, nplus>, Dbound);

System("mkdir -p tensor_indecomp");
Write("tensor_indecomp/" * sig, Join(indecomp, "\n"));
Write("tensor_indecomp/" * sig * ".done", "done");

System("mkdir -p lattice_decomp_data");
dname := "lattice_decomp_data/" * sig;
for label in indecomp do
    Write(dname, Sprintf("%o|\\N|f", label));
end for;
for label -> D in decomp do
    s := [];
    for opt in D do
        ctr := AssociativeArray(:Default:=0);
        for lab in opt do
            ctr[lab] +:= 1;
        end for;
        labs := [k : k->v in ctr];
        SortLabels(~labs);
        // labels quoted: the column is jsonb, and [[[2.2.1.5.1,1]]] is not JSON
        Append(~s, Join([Sprintf("[\"%o\",%o]", lab, ctr[lab]) : lab in labs], ","));
    end for;
    Dstr := "[[" * Join(s, "],[") * "]]";
    Write(dname, Sprintf("%o|%o|t", label, Dstr));
end for;
