intrinsic GetThetaBound(r::RngIntElt, N::RngIntElt, D::RngIntElt, is_even::BoolElt) -> RngIntElt
{Dimension of the modular-forms space containing the theta series of a rank-r,
level-N, determinant-D lattice with the given even/odd parity; used to bound how many
theta coefficients can distinguish lattices in a genus.}
    k := r div 2;
    Dstar := (-1)^k * D;
    kron := KroneckerCharacter(Dstar);
    // Theta of L lies in M_{r/2}(Gamma_0(base), chi): base = N when L is even of even
    // rank, else 4N (odd rank => half-integral weight; an odd lattice => level x 4).
    // Keying only on IsEven(r) is wrong for an odd lattice of even rank -- its character
    // can have conductor 4N > N, which then fails to coerce into the mod-N Dirichlet
    // group ("character has too large conductor").  The LCM with the actual conductor is
    // a belt-and-braces guard so no nonstandard conductor can re-trigger that crash.
    base := (is_even and IsEven(r)) select N else 4*N;
    X := DirichletGroupFull(LCM(base, Conductor(kron)));
    return Dimension(ModularForms(X!kron, k));
end intrinsic;

intrinsic HashGenus(G::SymGen) -> RngIntElt
{Compute the hash value associated to the genus of the given lattice}
    // Hash the canonical LMFDB genus label, not ConwaySymbol(G): equal genera can
    // have different Conway-symbol strings (e.g. a 2-adic oddity prints as "_1"
    // from Genus(L) but "_{I}" from GenusSymbolFromLabel), which broke the
    // label -> hash round-trip.  CreateGenusLabel is canonical, so
    // HashGenus(Genus(L)) = HashGenus(GenusSymbolFromLabel(LabelOf(L))).
    return CollapseIntList(StringToBytes(CreateGenusLabel(G)));
end intrinsic;

intrinsic HashGenus(L::Lat) -> RngIntElt
{Compute the hash value associated to the genus of the given lattice}
    return HashGenus(Genus(L));   // the custom SymGen hash, matching what is stored
end intrinsic;

intrinsic ThetaHash(theta_series::SeqEnum[RngIntElt], genus_hash::RngIntElt, prec::RngIntElt) -> RngIntElt
{Combines the hash of the genus with the truncation of the theta series with precision prec}
    // Note: should handle cases where theta_series has length less than prec (trailing zeros)
    // and more than prec (truncate)
    return CollapseIntList([genus_hash] cat theta_series[1..Min(prec,#theta_series)] cat [0 : _ in [1..prec-#theta_series]]);
end intrinsic;

intrinsic BVHash(L::Lat, genus_hash::RngIntElt, m::RngIntElt) -> RngIntElt
{Combines the hash of the genus with a BV hash of vectors of norm at most m}
    lat_hash := HBV(GramMatrix(L), m);
    return CollapseIntList([genus_hash, lat_hash]);
end intrinsic;

intrinsic BVhashes(lattices::SeqEnum[Lat], genus_hash::RngIntElt, m::RngIntElt) -> SeqEnum[RngIntElt]
{Multiple BVHashes at once, for use in TimeoutCall}
    return [BVHash(L, genus_hash, m) : L in lattices];
end intrinsic;

intrinsic HashLat(L::Lat, genus_hash::RngIntElt, hash_func::MonStgElt) -> RngIntElt
{Compute the hash value associated to the given lattice}
    code := hash_func[1..2];
    m := StringToInteger(hash_func[3..#hash_func]);
    if code eq "Th" then
        theta := Eltseq(ThetaSeries(L, m-1));
        return ThetaHash(theta, genus_hash, m);
    elif code eq "BV" then
        return BVHash(L, genus_hash, m);
    else
        error Sprintf("Invalid hash code %o", hash_func);
    end if;
end intrinsic;

function FirstDiff(f1, f2, M)
    // The smallest integer where f1 and f2 differ, or M if there is no difference up to precision M
    m := 0;
    while m lt M do
        c1 := (m lt #f1) select f1[m+1] else 0;
        c2 := (m lt #f2) select f2[m+1] else 0;
        if c1 ne c2 then
            return m;
        end if;
        m +:= 1;
    end while;
    return M;
end function;

function PickBest(hash_opts)
    // Given a sequence of triples <-distinguished, total_time, func_name>, pick the best
    // For now we just sort and take the first, but perhaps we should balance time and number of distinguished lattices rather than sorting lexicographically
    return Minimum(hash_opts);
end function;

intrinsic SetHashes(~lats::SeqEnum[Assoc], ~genus::Assoc, theta_elapsed::Assoc, timeout::RngIntElt)
{Set hash values and theta distinguishing information}
    L0 := lats[1]["lattice"];
    genus_hash := HashGenus(Genus(L0)); // TODO: could have problems if gram is too big to store in shortint[] so we're using gram_others.
    genus["genus_hash"] := genus_hash;

    // These are copied verbatim (as strings) from the genus's basic record, so
    // coerce to integers before any arithmetic use -- in particular GetThetaBound
    // below wants RngIntElt and otherwise crashes with "Bad argument types" (a bug
    // that only surfaces on the rare genus whose theta series do NOT distinguish it,
    // the one code path that reaches GetThetaBound).
    rank := StringToInteger(lats[1]["rank"]);
    nplus := StringToInteger(lats[1]["nplus"]);
    if rank ne nplus then
        // We don't compute hashes for indefinite lattices
        // TODO: maybe we should in rank 3, using spinor genera
        genus["theta_distinguishing_prec"] := "\\N";
        genus["is_theta_distinguished"] := "\\N";
        genus["hash_function"] := "\\N";
        genus["is_hash_distinguished"] := "\\N";
        for i in [1..#lats] do
            lats[i]["hash"] := "\\N";
        end for;
        return;
    end if;
    level := StringToInteger(lats[1]["level"]);
    det := StringToInteger(lats[1]["disc_abs"]);   // there is no "det" column; disc_abs is |determinant|

    thetas := Sort([<lat["theta_series"], lat["theta_prec"]> : lat in lats]);
    dprec := 0;
    dcount := 1;
    mprec := Minimum([lat["theta_prec"] : lat in lats]);
    mcount := 1;
    distinguished := true;
    for i in [1..#thetas-1] do
        M := Minimum(thetas[i][2], thetas[i+1][2]);
        m := FirstDiff(thetas[i][1], thetas[i+1][1], M);
        if m eq M then
            distinguished := false;
        else
            dprec := Maximum(dprec, m+1);
            dcount +:= 1;
            if m lt mprec then
                mcount +:= 1;
            end if;
        end if;
    end for;
    genus["theta_distinguishing_prec"] := dprec;
    if distinguished then
        genus["is_theta_distinguished"] := true;
    elif dprec ge GetThetaBound(rank, level, det, lats[1]["is_even"] eq "T") then
        genus["is_theta_distinguished"] := false;
    else
        genus["is_theta_distinguished"] := "\\N";
    end if;

    hash_opts := [];
    // Every "elapsed" that goes into a hash_opts tuple must land in ONE real field:
    // hash_opts has no declared universe, so its first Append fixes the tuple's
    // element types, and theta_elapsed mixes reals (compute times) with Infinity()
    // (type Infty, filled in for precisions a representative never reached).  If a
    // theta option with an Infinity elapsed is appended first, the 2nd slot becomes
    // Infty and the subsequent real BV elapsed can no longer coerce in -- the
    // "Append: Cannot coerce element into the universe" crash.  to_secs collapses
    // both cases to RealField (Infinity -> a large sentinel, since elapsed is only a
    // tie-breaker in PickBest = Minimum, where slower/unavailable should sort last).
    RR := RealField();
    to_secs := func<x | (x cmpeq Infinity()) select (RR ! 10^30) else (RR ! x)>;

    // If theta hash is an option, record how good it is
    if dprec le mprec then
        // Can use theta series for hash function
        tprec := Minimum([k : k in Keys(theta_elapsed) | k ge dprec]);
        Append(~hash_opts, <-dcount, to_secs(theta_elapsed[tprec]), Sprintf("Th%o", dprec)>);
    elif mprec gt 1 then
        tprec := Minimum([k : k in Keys(theta_elapsed) | k ge mprec]);
        Append(~hash_opts, <-mcount, to_secs(theta_elapsed[tprec]), Sprintf("Th%o", mprec)>);
    end if;

    // Now check BV options
    // "\N" where the minimum timed out (FillGenus computes it under a TimeoutCall)
    minima := [Integers() | lat["minimum"] : lat in lats | Type(lat["minimum"]) eq RngIntElt];
    min_norm := #minima gt 0 select Minimum(minima) else 3;
    if min_norm lt 3 then min_norm := 3; end if;
    BVvals := AssociativeArray();
    lattices := [lat["lattice"] : lat in lats];
    for m in [min_norm..min_norm+3] do
        // TODO: if there are too many vectors of norm up to m we can bail before calling BVhashes
        // HardKill: BVhashes enumerates vectors with ShortVectors, which Alarm cannot interrupt
        success, res, elapsed := TimeoutCall(timeout, BVhashes, <lattices, genus_hash, m>, 1 : HardKill := true);
        if success then
            res := res[1];   // TimeoutCall wraps the results in a List; unwrap to the hash sequence
            dcount := #{h : h in res};
            Append(~hash_opts, <-dcount, to_secs(StringToReal(elapsed)), Sprintf("BV%o", m)>);
            BVvals[m] := res;
            if dcount eq #lats then
                break;
            end if;
        else
            break;
        end if;
    end for;

    // Are there other invariants we should use?  root system, dual theta series...

    // TODO: Need to handle the case where all BV options timed out and some theta series completely failed (and thus hash_opts is empty)
    best := PickBest(hash_opts);
    genus["is_hash_distinguished"] := (best[1] eq -#lats);
    genus["hash_function"] := best[3];
    code := best[3][1..2];
    m := StringToInteger(best[3][3..#best[3]]);   // #best[3] (string length), not #best (tuple arity)
    if code eq "BV" then
        vals := BVvals[m];
    elif code eq "Th" then
        vals := [ThetaHash(lat["theta_series"], genus_hash, m) : lat in lats];
    else
        error "Invalid hash code";
    end if;

    for i in [1..#lats] do
        lats[i]["hash"] := vals[i];
    end for;
end intrinsic;


