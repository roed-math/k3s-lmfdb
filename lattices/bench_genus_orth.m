// Benchmark the orbit-method genus enumeration against the p-neighbour walk.
//
// Usage:
//   magma -b cases:=14.3:18.3:22.3 method:=both bench_genus_orth.m
//
// Options:
//   cases   colon-separated case names (default: a moderate ladder); see below
//   method  orth | faster | both (default both)
//   ftimeout  timeout in seconds for each GenusRepresentativesFaster run (default 3600)
//   verbose   verbosity for the GenusReps flag (default 1)
//
// Case names n.p refer to the genus G_{n,p} of even lattices of rank n and
// determinant p (n even) or 2p (n odd), seeded from Chenevier-Taibi Table 7.2
// padded with copies of E8; expected class numbers are their Tables 1.2/1.3.

AttachSpec("lattices.spec");

if not assigned verbose then verbose := "1"; end if;
SetVerbose("GenusReps", StringToInteger(verbose));
if not assigned method then method := "both"; end if;
if not assigned ftimeout then ftimeout := "3600"; end if;
ftimeout := StringToInteger(ftimeout);

Z := Integers();
CM := func<s | LatticeWithGram(CartanMatrix(s))>;
DS := function(lats)
    L := lats[1];
    for i in [2..#lats] do L := DirectSum(L, lats[i]); end for;
    return L;
end function;

// Even overlattice of index Order(x) obtained by adjoining an isotropic element
// of the discriminant group (used to build F8, the even det-5 rank-8 lattice).
function even_overlattice(L)
    LD := Dual(L : Rescale := false);
    A, quo := LD / L;
    n := Rank(L);
    for x in A do
        o := Order(x);
        if o le 1 then continue; end if;
        xl := x @@ quo;
        if not IsIntegral(Norm(xl)) or not IsEven(Z ! Norm(xl)) then continue; end if;
        M := VerticalJoin(ScalarMatrix(Rationals(), n, 1),
                          Matrix(Rationals(), 1, n, Eltseq(xl)));
        H := HermiteForm(ChangeRing(o * M, Z));
        B := ChangeRing(Matrix([H[i] : i in [1..n]]), Rationals()) / o;
        P := Lattice(B, ChangeRing(GramMatrix(L), Rationals()));
        GP := GramMatrix(P);
        if forall{ e : e in Eltseq(GP) | IsIntegral(e) } and IsEven(LatticeWithGram(ChangeRing(GP, Z))) then
            return LatticeWithGram(ChangeRing(GP, Z));
        end if;
    end for;
    error "no even overlattice found";
end function;

S2 := LatticeWithGram(Matrix(Z, 2, 2, [2,-1,-1,4]));
F8 := even_overlattice(DS([CM("E7"), LatticeWithGram(Matrix(Z,1,1,[10]))]));
assert Determinant(F8) eq 5 and IsEven(F8);

builders := AssociativeArray();
expected := AssociativeArray();
procedure addcase(~builders, ~expected, name, lats, h)
    builders[name] := lats;
    expected[name] := h;
end procedure;

E8 := CM("E8"); E6 := CM("E6"); A1 := CM("A1"); A2 := CM("A2");
A4 := CM("A4"); A5 := CM("A5"); A6 := CM("A6");

// p = 3, even rank (det 3)
addcase(~builders, ~expected, "14.3", [E6, E8], 2);
addcase(~builders, ~expected, "18.3", [A2, E8, E8], 6);
addcase(~builders, ~expected, "22.3", [E6, E8, E8], 31);
addcase(~builders, ~expected, "26.3", [A2, E8, E8, E8], 678);
// p = 3, odd rank (det 6)
addcase(~builders, ~expected, "11.3", [A1, A2, E8], 2);
addcase(~builders, ~expected, "15.3", [A1, E6, E8], 5);
// NB: the paper's Table 1.3 prints 19 for (19,3), but the authors' data site
// (https://olitb.net/pro/uni29/genera/even_hdiscp/3/X19.html) and our own
// mass-verified enumerations (orth + p-neighbours + built-in isometry tests)
// all give 27; the table entry is a typo.
addcase(~builders, ~expected, "19.3", [A1, A2, E8, E8], 27);
addcase(~builders, ~expected, "21.3", [A5, E8, E8], 64);
addcase(~builders, ~expected, "23.3", [A1, E6, E8, E8], 290);
// p = 5
addcase(~builders, ~expected, "12.5", [A4, E8], 2);
addcase(~builders, ~expected, "16.5", [F8, E8], 5);
addcase(~builders, ~expected, "20.5", [A4, E8, E8], 27);
addcase(~builders, ~expected, "24.5", [F8, E8, E8], 352);
addcase(~builders, ~expected, "13.5", [A1, A4, E8], 5);
addcase(~builders, ~expected, "21.5", [A1, A4, E8, E8], 210);
// p = 7
addcase(~builders, ~expected, "14.7", [A6, E8], 4);
addcase(~builders, ~expected, "18.7", [S2, E8, E8], 20);
addcase(~builders, ~expected, "22.7", [A6, E8, E8], 153);
addcase(~builders, ~expected, "15.7", [A1, A6, E8], 14);
addcase(~builders, ~expected, "19.7", [A1, S2, E8, E8], 119);
// composite determinant (det 27, disc group Z/3^3): chain 27 -> 9 -> 3 -> 1
addcase(~builders, ~expected, "14.27", [A2, A2, A2, E8], 0);
// odd lattice, det 3
addcase(~builders, ~expected, "I13A2", [LatticeWithGram(IdentityMatrix(Z, 13)), A2], 0);

if not assigned cases then
    cases := "11.3:14.3:12.5:13.5:15.3:16.5:14.7:15.7:18.3:18.7:19.3:14.27:I13A2";
end if;
case_list := Split(cases, ":");

function orth_wrap(L)
    return GenusRepresentativesOrth(L);
end function;
function faster_wrap(L)
    return GenusRepresentativesFaster(L);
end function;

printf "case       rank det  even    h  orth(s)  orth2(s)  faster(s)\n";
for name in case_list do
    error if not IsDefined(builders, name), "unknown case " cat name;
    L := DS(builders[name]);
    h := expected[name];
    n := Rank(L); d := Determinant(L);
    horth := -1; torth := -1.0; torth2 := -1.0;
    hfast := -1; tfast := -1.0;
    if method in ["orth", "both"] then
        t0 := Cputime();
        reps := GenusRepresentativesOrth(L);
        torth := Cputime(t0);
        horth := #reps;
        error if h ne 0 and horth ne h,
            Sprintf("case %o: orth found %o classes, expected %o", name, horth, h);
        error if not forall{ R : R in reps | Genus(R) eq Genus(L) },
            Sprintf("case %o: orth produced a lattice outside the genus", name);
        // warm-cache rerun (top-level cache hit; shows cache overhead only)
        t0 := Cputime();
        _ := GenusRepresentativesOrth(L);
        torth2 := Cputime(t0);
    end if;
    if method in ["faster", "both"] then
        okf, res, elapsed := TimeoutCall(ftimeout, faster_wrap, <L>, 1);
        if okf then
            hfast := #res[1];
            tfast := StringToReal(Sprintf("%o", elapsed));  // TimeoutCall reports elapsed as a string
            error if h ne 0 and hfast ne h,
                Sprintf("case %o: faster found %o classes, expected %o", name, hfast, h);
        else
            printf "case %o: faster timed out after %o s\n", name, ftimeout;
        end if;
    end if;
    error if horth ge 0 and hfast ge 0 and horth ne hfast,
        Sprintf("case %o: orth found %o classes but faster found %o", name, horth, hfast);
    R3 := RealField(4);
    printf "%o %o %o %o %o %o %o %o\n",
        name, n, d, IsEven(L) select "even" else "odd", Maximum(horth, hfast),
        R3 ! torth, R3 ! torth2, R3 ! tfast;
end for;
printf "peak memory: %o MB\n", GetMaximumMemoryUsage() div 2^20;
exit 0;
