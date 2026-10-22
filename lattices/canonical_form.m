// freeze;

declare verbose CanonicalForm, 5;

function IsWellRounded(A)
    L := LatticeWithGram(A);
    minA := ShortestVectors(L);
    LminA := Lattice(ChangeRing(Matrix(minA),Rationals()), A);
    return Rank(LminA) eq Rank(L);
end function;

function V_wr_cvp(A : bound := Infinity())
    assert IsWellRounded(A);
    L := LatticeWithGram(A);
    minA := &cat[[v,-v] : v in ShortestVectors(L)];
    LminA := Lattice(ChangeRing(Matrix(minA),Rationals()), A);
    assert Rank(LminA) eq Rank(L);
    Zn := LatticeWithGram(A);
    Q, quo := Zn / LminA;
    cs := [ChangeRing(Vector(q@@quo), Rationals()) : q in Q];
    B := ChangeRing(Matrix(Basis(LminA)),Rationals());
    L_B := LatticeWithGram(B * A * Transpose(B));
    ret := [Vector(ChangeRing(v, Rationals())) : v in minA];
    if Type(bound) eq RngIntElt then
        ret cat:= &cat[[c - Vector(ChangeRing(v, Rationals())*B) : v in ClosestVectors(L_B, c*B^(-1) : Max := bound)] : c in cs | not IsZero(c)];
    else
        ret cat:= &cat[[c - Vector(ChangeRing(v, Rationals())*B) : v in ClosestVectors(L_B, c*B^(-1))] : c in cs | not IsZero(c)];
    end if;
    return ret;
end function;

function satspan(v, A)
    v_mat := Matrix(v);
    d := Denominator(v_mat);
    dv_mat := ChangeRing(d*v_mat, Integers());
    S, X, Y := SmithForm(dv_mat);
    one := IdentityMatrix(Integers(), Nrows(S));
    zero := ZeroMatrix(Integers(), Nrows(S), Ncols(S)-Nrows(S));
    B := HorizontalJoin(one, zero);
    sat_basis := ChangeRing(X^(-1)*B*Y^(-1), Rationals());
    return Lattice(sat_basis, A);
end function;

function V_ms(A : max_num := Infinity())
    vprintf CanonicalForm, 4:
            "\n\t\t\t\t Creating characteristic set from short vectors...";
    A := ChangeRing(A, Rationals());
    L := LatticeWithGram(A : CheckPositive := false);
    max_norm := AbsoluteValue(Maximum(Diagonal(A)));
    // This takes too long, and we might not need all of them
    // all_vs := ShortVectors(L, max_norm);
    // norms := {v[2] : v in all_vs};
    // norms := Sort([n : n in norms]); 
    VA := [];
    norms := IsEven(L) select [2..max_norm by 2] else [1..max_norm];
    for n in norms do
        vprintf CanonicalForm, 5: "\n\t\t\t\t\t adding vectors of length %o...", n;
        if Type(max_num) eq RngIntElt then
            VA cat:= &cat[[v[1],-v[1]] : v in ShortVectors(L,n,n : Max := (max_num - #VA) div 2)];
        else
            VA cat:= &cat[[v[1],-v[1]] : v in ShortVectors(L,n,n)];
        end if;
        Lsub := sub<L | VA>;
        vprintf CanonicalForm, 5: "done! rank is %o.", Rank(Lsub);
        if Rank(Lsub) eq Rank(L) then
            vprintf CanonicalForm, 5: "quotient is %o.", L / Lsub;
        end if;
        if (Lsub eq L) then return VA; end if;
        if (#VA ge max_num) then return []; end if;
    end for;
    Lsub := sub<L | VA>;
    assert Lsub eq L;
    vprintf CanonicalForm, 4 : "done!\n";
    return VA;
end function;

function V_cvp(A : max_num := Infinity(), bound := Infinity())
    A := ChangeRing(A, Rationals());
    L := LatticeWithGram(A);
    if Type(max_num) eq RngIntElt then
        minA := &cat[[v,-v] : v in ShortestVectors(L : Max := max_num div 2)];
    else
        minA := &cat[[v,-v] : v in ShortestVectors(L)];
    end if;
    LminA := sub<L | minA>;
    if (LminA eq L) then
	    return minA;
    end if;
    if (#minA ge max_num) then
       return [];
    end if;
    L1 := satspan(Basis(LminA), A);
    assert LminA subset L1;
    assert IsFree(L/L1); // verifying L1 is saturated
    B1 := ChangeRing(Matrix(Basis(L1)),Rationals());
    A1 := B1 * A * Transpose(B1);
    B2 := Transpose(Matrix(Basis(Kernel(A*Transpose(B1)))));
    proj := B2 * (Transpose(B2)*A*B2)^(-1) * Transpose(A*B2);
    B2 := ChangeRing(Matrix(Basis(Lattice(Transpose(proj)))), Rationals());
    A2 := B2 * A * Transpose(B2);
    r := Rank(L1);
    if (r eq Rank(L)) then
	    V_cvp_A2 := [];
    else
	    V_cvp_A2 := V_cvp(A2 : max_num := max_num, bound := bound);
        if IsEmpty(V_cvp_A2) then return []; end if;
    end if;
    A1_part := [Vector(Rationals(), v)*B1 : v in V_wr_cvp(A1 : bound := bound)];
    proj_Z := ChangeRing(Denominator(proj)*proj, Integers());
    assert sub<L | A1_part> eq L1;
    // This is not true, need to figure out what the correct statement is
    // assert RSpaceWithBasis(BasisMatrix(L1)) eq Kernel(proj_Z);
    B2_Z := ChangeRing(Denominator(proj)*B2, Integers());
    A2_part := [Solution(Transpose(proj_Z), Vector(v)*B2_Z) : v in V_cvp_A2];
    // union_A2_part := &cat[[v - w : w in ClosestVectors(L1,v - v*proj)] : v in A2_part];
    if Type(bound) eq RngIntElt then
        union_A2_part :=  &cat[[v - w : w in ClosestVectors(L1,v - v*Transpose(proj) : Max := bound)] : v in A2_part];
    else
        union_A2_part :=  &cat[[v - w : w in ClosestVectors(L1,v - v*Transpose(proj))] : v in A2_part];
    end if;
    VA := A1_part cat union_A2_part;
    VA := [v : v in VA | not IsZero(v)];
    Lsub :=  sub<L | VA>;
    if Type(bound) ne RngIntElt then
        assert Lsub eq L;
    end if;
    return VA;
end function;

function V_best(A)
    // Ad := GramMatrix(LLL(Dual(LatticeWithGram(A))));
    // Ad := GramMatrix(Dual(LatticeWithGram(A) : Rescale := false));
    VAs := [];
    // At the moment something is wrong with using duals, so we don't use them
    // for B in [A, Ad] do
    for B in [A] do
	    Append(~VAs, <[Vector(v) : v in V_ms(B)],B>);
        try
            Append(~VAs, <[Vector(v) : v in V_cvp(B )],B>);
        catch e
            continue;
        end try;
    end for;
    sorted := Sort(VAs, func<x,y | #x[1]-#y[1]>);
    // return sorted[1];
    return sorted[1][1];
end function;

function minimize_vector_set(VA, G, L)
    //G := AutomorphismGroup(Dual(badl));
    //G := AutomorphismGroup(L);
    lut := AssociativeArray();
    for j->vj in VA do lut[vj] := j; end for;
    perms := [[lut[VA[j]^g] : j in [1..#VA]] : g in Generators(G)];
    pg := PermutationGroup<#VA|perms>;
    orbs := Orbits(pg);
    pi := [#orb : orb in orbs];
    subs := [S : S in Subsets({1..#pi})];
    lens := [&+[Integers() | pi[j] : j in S] : S in subs];
    Sort(~lens, ~orb_subs_perm);
    for j in [1..#subs] do
        S := subs[j^orb_subs_perm];
        vecs := &cat[[VA[i] : i in orbs[j]] : j in S];
        sublat := sub< L | vecs>;
        if Rank(sublat) ne Rank(L) then continue; end if;
        if sublat eq L then return vecs; end if;
        // Not yet sure what to do if we only span a sublattice of finite index
        /*
        A := ChangeRing(GramMatrix(L), Rationals());
        L1 := satspan(Basis(sublat), A);
        if L1 eq L then
            B1 := ChangeRing(Matrix(Basis(L1)),Rationals());
            A1 := B1 * A * Transpose(B1);
            vecs := [Vector(Rationals(), v)*B1 : v in V_cvp(A1)];
            assert sub<L | vecs> eq L;
            return vecs;
        end if;
        */
    end for;
    // Should not be reaching this point
    assert false;
end function;

// Largest characteristic vector set CanonicalForm will work with (see V_best_with_dual).
MaxCharVectors := 4096;

function V_best_with_dual(A)
    // We want to track the transformations, so we create the dual by hand
    // Ld, Td := LLL(Dual(LatticeWithGram(A)));
    // Ad := GramMatrix(Ld);
    Ainv := Determinant(A)*A^(-1);
    Ld, Td := LLL(LatticeWithGram(Ainv));
    Ad := GramMatrix(Ld);
    assert Ad eq Td*Ainv*Transpose(Td);
    mats := [A, Ad];
    is_dual := [false, true];
    VAs := [];
    // Cap on the characteristic vector set.  Everything downstream is quadratic in
    // its size (W = B*A*B^t and T1's matrix are #VA x #VA over Q, then a canonical
    // graph labelling on ~#VA*log #VA vertices), and the cap used to start at 2^20
    // and double: on the C=32768 census a rank-17 lattice reached 33250 vectors, a
    // 1.1e9-entry matrix, with Magma processes at 130-220 GB.  Past the cap we fail
    // at once -- the same outcome as a CanonicalForm timeout (no canonical gram) --
    // instead of allocating it.
    max := MaxCharVectors;
    found := false;
    vprintf CanonicalForm, 1: "Finding best characteristic vector set...";
    while not found do
        vprintf CanonicalForm, 2: "\n\t allowing for up to %o vectors", max;
        for j->B in mats do
            vprintf CanonicalForm, 3: "\n\t\t trying the %olattice...", is_dual[j] select "dual " else "";
            for jj->Vchar in [V_ms, V_cvp] do
                vprintf CanonicalForm, 4: "\n\t\t\t trying the %o characteristic vector set...", jj eq 1 select "V_ms" else "V_cvp";
		        VA := [Vector(v) : v in Vchar(B : max_num := max)];
                if not IsEmpty(VA) then
                    Append(~VAs, <VA,B,is_dual[j], Td>);
                    if (#VA lt max) then
		                found := true;
                        max := #VA;
                    end if;
                end if;
                vprintf CanonicalForm, 4 : "Done!"; 
            end for;
            vprintf CanonicalForm, 3 : "Done!";
         end for;
         vprintf CanonicalForm, 2 :"Done!"; 
         error if not found, Sprintf("CanonicalForm: no characteristic vector set of at most %o vectors", MaxCharVectors);
    end while;
    vprintf CanonicalForm, 1 : "Done!\n";
    sorted, perm := Sort(VAs, func<x,y | #x[1]-#y[1]>);
    vprintf CanonicalForm, 3 : "Choosing method %o\n", 1^perm;
    return sorted[1];
end function;

function T1(W)
    a := Maximum(Eltseq(W))+1;
    b := a + 1;
    p := Nrows(W);
//    W_prime := MatrixRing(Integers(),p+2)!0;
    W_prime := MatrixRing(Rationals(),p+2)!0;
    for i in [1..p] do
	for j in [i+1..p] do
	    W_prime[i,j] := W[i,j];
	end for;
	W_prime[i,p+1] := W[i,i];
	W_prime[i,p+2] := a;
    end for;
    W_prime[p+1,p+2] := b;
    for i in [1..p+2] do
	for j in [1..i-1] do
	    W_prime[i,j] := W_prime[j,i];
	end for;
    end for;
    return W_prime;
end function;

function bits(s,w)
    return [BitwiseAnd(ShiftRight(s,k),1) : k in [0..w-1]];
end function;

function T2(W)
    // Apparently the sort key doesn't really matter,
    // as the labels will be determined by the binary digits.
    // S := Sort([s : s in Set(Eltseq(W))], 
    // func< x,y | Abs(x) eq Abs(y) select x-y else Abs(x)-Abs(y) >);
    S := Sort([s : s in Set(Eltseq(W))]);
    w := Ceiling(Log(2,#S));
    all_bits := [bits(s,w) : s in [0..#S-1]];
    all_bits_bool := [[all_bits[s][k] eq 1 : k in [1..w]]
		      : s in [1..#S]];
    p := Nrows(W);
    lut := AssociativeArray(S);
    for i->s in S do
	lut[s] := i;
    end for;
    vertices := [];
    edges := [];
    for i in [1..p] do
	v_nbrs := [[<i_prime,k-1> : i_prime in [i+1..p] | all_bits_bool[lut[W[i,i_prime]]][k] ] cat [<i,k_prime> : k_prime in [k..w-1]] : k in [1..w]];
	vertices cat:= [<i,k> : k in [0..w-1]];
	edges cat:= v_nbrs;
    end for;
    vertices := {@ v : v in vertices @};
    edges := [{@ e : e in es @} : es in edges];
    wts := [v[2] : v in vertices];
    G := Graph<vertices|edges>;
    AssignLabels(VertexSet(G),wts);
    return G;
end function;

function U_V(A)
    A := ChangeRing(A, Rationals());
    VA, A, is_dual, Td := Explode(V_best_with_dual(A));
    vprintf CanonicalForm, 4 : "VA is %o\n", VA;
    // VA := V_best(A);
    // This is not really needed, we just keep track of the weights
    // G_A := CompleteGraph(#VA);
    B := ChangeRing(Matrix(VA),Rationals());
    W := B * A * Transpose(B);
    G := T1(W);
    T2G := T2(G);
    T2G_can, _, _, labels := CanonicalGraphMap(T2G);
    p := Nrows(G);
    S := [s : s in Set(Eltseq(G))];
    w := Ceiling(Log(2,#S));
    perm := [Minimum([Index(labels, <i,k>) : k in [0..w-1]]) : i in [1..p]];
    // this is the ordering of the vertices in G
    perm_inv := [Index(perm, i) : i in [1..#perm]];
    p -:= 2;
    v := [VA[i] : i in perm_inv | i le p];
    QA := ChangeRing(Transpose(Matrix(v)), Integers());
    H, U_inv := HermiteForm(QA);
    U := U_inv^(-1);
    assert U*H eq Transpose(Matrix(v));
    vprintf CanonicalForm, 4 : "U is %o\n", U;
    return U, A, is_dual, Td;
end function;

/*
function CanonicalForm(A)
    U := U_V(A);
    can_A := Transpose(U)*A*U;
    return can_A;
end function;
*/

intrinsic CanonicalForm(A::AlgMatElt) -> AlgMatElt
{Return a canonical form for the definite matrix A, together with a matrix T such that Acan = T*A*Transpose(T)}
    if IsNegativeDefinite(A) then
      mAcan, mU := CanonicalForm(-A);
      return -mAcan, mU;
    end if;
    require IsPositiveDefinite(A) : "A must be definite";
    vprintf CanonicalForm, 1 : "Computing canonical form for \n %o\n", A;
    L, T := LLL(LatticeWithGram(A : CheckPositive := false));
    Ared := GramMatrix(L);
    assert T*A*Transpose(T) eq Ared;
    U, Ad, is_dual, Td := U_V(Ared);
    can_A := Transpose(U)*Ad*U;
    if is_dual then
      Ared := ChangeRing(Ared, Rationals());
      A := ChangeRing(A, Rationals());
      assert Ad eq Td*Determinant(Ared)*Ared^(-1)*Transpose(Td);
      U1 := Transpose(U)*Td;
      Aredd := Determinant(Ared)*Ared^(-1);
      assert can_A eq U1*Aredd*Transpose(U1);
      assert Aredd eq Determinant(T)^2*Determinant(A)*Transpose(T)^(-1)*A^(-1)*T^(-1);
      U2 := U1*Transpose(T)^(-1);
      D := Determinant(T)^2*Determinant(A);
      Ainv := D*A^(-1);
      assert can_A eq U2*Ainv*Transpose(U2);
      U3 := Transpose(U2)^(-1);
      assert can_A^(-1) eq U3*Ainv^(-1)*Transpose(U3);
      assert D*can_A^(-1) eq U3*A*Transpose(U3);
      can_A := D*can_A^(-1);
      U := U3;
    else
      U := Transpose(U)*T;
    end if;
    vprintf CanonicalForm, 4 : "can_A before sorting is %o\n", can_A;
    diag := Diagonal(can_A);
    Sort(~diag, ~perm);
    P := PermutationMatrix(Integers(), perm);
    U := P*U;
    can_A := P*can_A*Transpose(P);
    /*
    vprintf CanonicalForm, 4 : "can_A after sorting, before last LLL is %o\n", can_A;
    assert can_A eq U*A*Transpose(U);
    L, T := LLL(LatticeWithGram(can_A : CheckPositive := false));
    can_A := GramMatrix(L);
    U := T*U;
    assert U*A*Transpose(U) eq can_A;
    */
    vprintf CanonicalForm, 1 : "Done!\n";
    vprintf CanonicalForm, 2 : "Canonical form is %o\n", can_A;
    return can_A, U;
end intrinsic;

intrinsic CanonicalForm(L::Lat) -> AlgMatElt
{Return a canonical form for L.}
  b := Basis(L);
  gram := InnerProductMatrix(L);
  b := [ChangeRing(Vector(Eltseq(x)), Rationals()) : x in b];
  A := Matrix(b)*gram*Transpose(Matrix(b));
  return CanonicalForm(A);
end intrinsic;

/*
intrinsic CanonicalForm(L::ModDedLat) -> AlgMatElt
{Return a canonical form for L.}
  lat := ZLattice(L);
  return CanonicalForm(lat);
end intrinsic;

intrinsic CanonicalForm(L::LatNF) -> AlgMatElt
{Return a canonical form for L.}
  return CanonicalForm(LatticeFromLatNF(L));
end intrinsic;
*/

// checks that V transforms well under linear transformations
procedure test_V(A : Vchar := V_cvp)
    A := ChangeRing(A, Rationals());
    T := RandomGLnZ(Nrows(A),33,17);
    T := ChangeRing(T, Rationals());
    TA := T*A*Transpose(T);
    A := GramMatrix(LLL(LatticeWithGram(A)));
    TA := GramMatrix(LLL(LatticeWithGram(TA)));
    dA := Denominator(A);
    _, T := IsIsometric(LatticeWithGram(dA*A), LatticeWithGram(dA*TA));
    T := ChangeRing(T, Rationals());
    assert TA eq T*A*Transpose(T);
    VA := Vchar(A);
    VTA := Vchar(TA);
    VA := [ChangeRing(v,Rationals()) : v in VA];
    VTA := [ChangeRing(v,Rationals()) : v in VTA];
    assert {w*T : w in VTA} eq Set(VA);
    return;
end procedure;

procedure test_V_with_dual(A : Vchar := V_best_with_dual)
    A := ChangeRing(A, Integers());
    T := RandomGLnZ(Nrows(A),33,17);
    TA := T*A*Transpose(T);
    A := GramMatrix(LLL(LatticeWithGram(A)));
    TA := GramMatrix(LLL(LatticeWithGram(TA)));
    _, T := IsIsometric(LatticeWithGram(A), LatticeWithGram(TA));
    assert TA eq T*A*Transpose(T);
    A := ChangeRing(A, Rationals());
    T := ChangeRing(T, Rationals());
    TA := ChangeRing(TA, Rationals());
    VB, B := Explode(Vchar(A));
    VTB, TB := Explode(Vchar(TA));
    _, T := IsIsometric(LatticeWithGram(B), LatticeWithGram(TB));
    T := ChangeRing(T, Rationals());
    VB := [ChangeRing(v,Rationals()) : v in VB];
    VTB := [ChangeRing(v,Rationals()) : v in VTB];
    assert {w*T : w in VTB} eq Set(VB);
    return;
end procedure;


// check that Stab(A) maps into Stab(G_A)
// Cannot establish surjectivity yet.
procedure test_W(A)
    A := ChangeRing(A, Integers());
    A := GramMatrix(LLL(LatticeWithGram(A)));
    A := ChangeRing(A, Rationals());
    VA := V_cvp(A);
    B := ChangeRing(Matrix(VA),Rationals());
    W := B * A * Transpose(B);
    stabA := AutomorphismGroup(LatticeWithGram(A));
    gens := Generators(stabA);
    // checking property (ii) of VA
    assert &and[{v*ChangeRing(S,Rationals()) : v in VA} eq Set(VA) : S in gens];
    imgs := [];
    for S in gens do
	sig := [Index(VA, v*ChangeRing(S, Rationals())) : v in VA];
	Append(~imgs, sig);
    end for;
    p := #VA;
    assert &and[&and[W[sig[i],sig[j]] eq W[i,j] : i,j in [1..p]] : sig in imgs];
end procedure;

function autT1_to_aut_T2(perm, w, vertexGset)
    p := #GSet(Parent(perm));
    big_perm := &cat[[Index(vertexGset, vertexGset[(i^perm-1)*w+j])
		      : j in [1..w]] : i in [1..p]];
    return Sym(w*p)!big_perm;
end function;

function autT2_to_autT1(sig, w, p)
    perm := [((w*i+1)^sig) div w + 1 : i in [0..p+1]];
    return Sym(p+2)!perm;
end function;

function autG_to_autT1(sig)
    p := #GSet(Parent(sig));
    perm := [i^sig : i in [1..p]] cat [p+1,p+2];
    return Sym(p+2)!perm;
end function;

function autT1_to_autG(sig)
    p := #GSet(Parent(sig))-2;
    perm := [i^sig : i in [1..p]];
    return Sym(p)!perm;
end function;
			 

procedure test_T2(A)
    A := ChangeRing(A, Integers());
    A := GramMatrix(LLL(LatticeWithGram(A)));
    A := ChangeRing(A, Rationals());
    VA := V_cvp(A);
    B := ChangeRing(Matrix(VA),Rationals());
    W := B * A * Transpose(B);
    G := T1(W);
    T2G := T2(G);
    autT2G, vertexGset, _, _, _, _ := AutomorphismGroup(T2G);
    act_vertices := Action(autT2G, vertexGset);
    p := #VA;
    S := [s : s in Set(Eltseq(G))];
    w := Ceiling(Log(2,#S));
    perms := [[Index(vertexGset, vertexGset[i]^act_vertices(sig))
	       : i in [1..(p+2)*w]] : sig in Generators(autT2G)];
    stabT2G := sub<Sym((p+2)*w) | perms>;
    stabA := AutomorphismGroup(LatticeWithGram(A));
    gens := Generators(stabA);
    imgs := [];
    for g in gens do
	assert {v*ChangeRing(g, Rationals()) : v in VA} eq Set(VA);
	sig := [Index(VA, v*ChangeRing(g, Rationals())) : v in VA];
	Append(~imgs, sig);
    end for;
    // imgs are now automorphisms of W, turning to automorphisms of G
    imgs := [sig cat [p+1,p+2] : sig in imgs];
    assert &and[&and[G[sig[i],sig[j]] eq G[i,j] : i,j in [1..p+2]]
		: sig in imgs];
    // converting to automorphisms of T2G
    big_perms := [&cat[[Index(vertexGset, vertexGset[(sig[i]-1)*w+j])
			: j in [1..w]] : i in [1..p+2]] : sig in imgs];
    stabA_im := sub<Sym((p+2)*w) | big_perms>;
end procedure;

function num_V(A)
    VA, A := Explode(V_best_with_dual(A));
    return #VA;
end function;

procedure test_canonical(A)
    T := Matrix(RandomGLnZ(Nrows(A),33,17));
    TA := T*A*Transpose(T);
    canA := CanonicalForm(A);
    canTA := CanonicalForm(TA);
    assert canA eq canTA;
end procedure;
