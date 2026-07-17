declare verbose FillGenus, 1;

import "neighbours.mag" : neighbours;

function hecke_primes(rank)
    if rank lt 8 then
        return [2,3,5];
    else
        return [2];
    end if;
end function;

intrinsic StringToReal(s::MonStgElt) -> RngIntElt
{ Converts a decimal string (like 123.456 or 1.23456e40 or 1.23456e-10) to a real number at default precision. }
    if #s eq 0 then return 0.0; end if;
    if "e" in s then
        t := Split(s,"e");
        require #t eq 2: "Input should have the form 123.456e20 or 1.23456e-10";
        return StringToReal(t[1])*10.0^StringToInteger(t[2]);
    end if;
    t := Split(s,".");
    require #t le 2: "Input should have the form 123 or 123.456 or 1.23456e-10";
    n := StringToInteger(t[1]);  s := t[1][1] eq "-" select -1 else 1;
    return #t eq 1 select RealField()!n else RealField()!n + s*RealField()!StringToInteger(t[2])/10^#t[2];
end intrinsic;

intrinsic ThetaSeriesIncremental(L::Lat, target_prec::RngIntElt, timeout::RngIntElt) -> SeqEnum, RngIntElt, Assoc
{Compute the theta series coefficients of L, doubling the precision until either target_prec is reached or the time budget timeout (in seconds) is exhausted. Returns the coefficient sequence obtained, the precision actually reached, and an associative array mapping each precision computed to the time (in seconds) it took to compute the theta series to that precision.}
    best_theta := [];
    best_prec := 0;
    theta_elapsed := AssociativeArray();
    remaining := timeout;
    prec := Maximum(16, Minimum(L) + 4);
    while prec le target_prec and remaining gt 0 do
        current_prec := Minimum(prec, target_prec);
        success, theta, elapsed := TimeoutCall(remaining, ThetaSeries, <L, current_prec - 1>, 1);
        if not success then
            vprintf FillGenus, 1 : "Theta series timed out at precision %o\n", current_prec;
            break;
        end if;
        best_theta := Eltseq(theta[1]);
        best_prec := current_prec;
        // ThetaSeries recomputes the whole series each time, so this elapsed time
        // is the cost of computing theta to current_prec (which is what SetHashes
        // uses to cost theta-based hashing).
        theta_elapsed[current_prec] := StringToReal(elapsed);
        vprintf FillGenus, 1 : "Theta series to precision %o in %o s\n", current_prec, elapsed;
        if current_prec ge target_prec then break; end if;
        remaining -:= Ceiling(StringToReal(elapsed));
        prec *:= 2;
    end while;
    return best_theta, best_prec, theta_elapsed;
end intrinsic;

function dict_to_jsonb(dict)
    return "{" * Join([Sprintf("\"%o\":%o", key, dict[key]) : key in Keys(dict)], ",") * "}";
end function;

intrinsic to_postgres(val::Any : jsonb_val := false) -> Any
{Serialise a value (matrix, sequence, tuple, associative array or scalar) into the Postgres array / JSON text format used for the database export.}
    delims := jsonb_val select "[]" else "{}";
    if ISA(Type(val),Mtrx) then
        return to_postgres(Eltseq(val) : jsonb_val:=jsonb_val);
    elif val cmpeq "\\N" then
        return val;
    // I think this is unnecessary, and used to cause problems, so removing for now.
    //elif Type(val) eq MonStgElt then
    //    return "\"" * val * "\""; // This will fail if the string has quotes, but I don't think that's ever true for us.
    elif Type(val) in [SeqEnum, Tup] then
        return delims[1] * Join([Sprintf("%o",to_postgres(x : jsonb_val:=jsonb_val)) : x in val],",") * delims[2];
    elif Type(val) eq Assoc then
        val_prime := AssociativeArray();
        for key in Keys(val) do
            val_prime[to_postgres(key)] := to_postgres(val[key] : jsonb_val:=true);
        end for;
        return dict_to_jsonb(val_prime);
    else
        return val;
    end if;
end intrinsic;

function RescaledDualNF(L)
    Q := Rationals();
    K := BaseRing(L);
    B := ChangeRing(BasisMatrix(L),Q);
    M := ChangeRing(InnerProductMatrix(L),Q);
    F := ChangeRing(GramMatrix(L),Q);
    B := F^-1 * B;
    B := IntegralMatrix(B);
    B div:= GCD(Eltseq(B));
    F := B * M * Transpose(B);
    F, d := IntegralMatrix(F);
    g := GCD(Eltseq(F));
    F div:= g;
    M := (d/g) * M;
    return NumberFieldLattice(Rows(ChangeRing(B, K)) : InnerProduct := ChangeRing(M,K));
end function;

function genus_reps_Magma(L)
    // The bound is set to infinity to avoid Magma printing an error message
    // without throwing a runtime error.
    if IsPositiveDefinite(GramMatrix(L)) or (Rank(L) eq 2) then
      return GenusRepresentatives(L : Bound := Infinity());
    end if;
    // due to some bugs in Magma, we convert to number field
    LF := NumberFieldLattice(L);
    reps := GenusRepresentatives(LF);
    return [LatticeWithGram(ChangeRing(GramMatrix(r), Integers()) :
			    CheckPositive := false) : r in reps];
end function;

function genus_reps_Logan(L)
    return Setseq(neighbours(L : thorough));
end function;

function genus_reps_Faster(L)
    // Our own mass-verified p-neighbour enumeration (definite lattices only).
    reps := GenusRepresentativesFaster(L);
    return reps;
end function;

function genus_reps_Spinor(L)
    // Indefinite lattices of rank >= 3: by Eichler each proper spinor genus is a single
    // class, so the spinor-genus representatives ARE the genus representatives.  This is
    // a fallback for when Magma's GenusRepresentatives fails (e.g. signature (2,2) rank 4
    // -- "Illegal null sequence" in its number-field-lattice code); SpinorRepresentatives
    // uses a different path and succeeds where GenusRepresentatives does not.
    reps := SpinorRepresentatives(L);
    return reps;
end function;

function SphereVolume(n)
    RR := RealField();
    pi := Pi(RR);
    m := n div 2;
    if IsEven(n) then
        return pi^m / Factorial(m);
    else
        return 2^n * pi^m * Factorial(m) / Factorial(n);
    end if;
end function;

// Exact GL2(Z)-isometry test for the canonical rank-2 forms [[0,m],[m,k1]] and
// [[0,m],[m,k2]] of (square) determinant -m^2.  These isotropic forms have finite
// isometry groups, so we can decide isometry directly: an isometry sends a
// primitive isotropic vector of the first to one of the second, and the rest of
// the basis is then determined up to the (finitely many) sign/line choices.
function square_disc_isometric(k1, k2, m)
    G1 := Matrix(Integers(), 2, 2, [0, m, m, k1]);
    G2 := Matrix(Integers(), 2, 2, [0, m, m, k2]);
    g1 := GCD(k1, 2*m);                       // GCD(0, 2m) = 2m
    isotropic := [ Vector(Integers(), [1, 0]),
                   Vector(Integers(), [-k1 div g1, (2*m) div g1]) ];
    for v1 in isotropic, sgn in [1, -1] do
        w1 := sgn * v1;
        row := Vector(Integers(), Eltseq(Matrix(w1) * G1));   // w1^t G1
        d := GCD(row[1], row[2]);
        if d eq 0 or m mod d ne 0 then continue; end if;
        _, p, q := XGCD(row[1], row[2]);                      // p*row[1] + q*row[2] = d
        v20 := (m div d) * Vector(Integers(), [p, q]);        // <w1, v20> = m
        Qv20 := (Matrix(v20) * G1 * Transpose(Matrix(v20)))[1][1];
        if (k2 - Qv20) mod (2*m) ne 0 then continue; end if;  // Q(v20 + t w1) = Qv20 + 2 t m
        v2 := v20 + ((k2 - Qv20) div (2*m)) * w1;
        U := Matrix(Integers(), 2, 2, [w1[1], v2[1], w1[2], v2[2]]);   // columns w1, v2
        if Abs(Determinant(U)) eq 1 and Transpose(U)*G1*U eq G2 then
            return true;
        end if;
    end for;
    return false;
end function;

// Genus representatives of a rank-2 lattice L0 of square determinant -m^2.  Magma's
// GenusRepresentatives fails here (it reduces to binary quadratic forms, which
// reject square discriminants).  Every such lattice is isometric to a canonical
// [[0,m],[m,k]] with 0 <= k < 2m; we keep those in L0's genus (Genus() is
// reliable for these even though the isometry routines are not) and remove
// duplicate isometry classes with the exact test above.
function genus_reps_square_disc(L0)
    m := Isqrt(Integers() ! (-Determinant(L0)));
    G0 := Genus(L0);
    reps := [];
    rep_ks := [];
    for k in [0 .. 2*m - 1] do
        Lk := LatticeWithGram(Matrix(Rationals(), 2, 2, [0, m, m, k]) : CheckPositive := false);
        if Genus(Lk) ne G0 then continue; end if;
        if forall{ kr : kr in rep_ks | not square_disc_isometric(k, kr, m) } then
            Append(~rep_ks, k);
            Append(~reps, Lk);
        end if;
    end for;
    return reps;
end function;

intrinsic FillGenus(label::MonStgElt : timeout := 1800, masslimit := 0, sizelimit := 0, timelimit := 0, adjlimit := 0)
{Fill the data for a genus and its lattice representatives, given files in the genera_basic format.

Enumeration guards (0 = unlimited): masslimit skips enumerating a definite genus whose
mass exceeds it (its class number is too large to enumerate); adjlimit skips the
adjacency (Hecke) matrix once the class number exceeds it (its cost is ~class_number^2);
sizelimit records the genus-level data but skips storing individual lattices past that
class number; timelimit is a per-genus wall-clock cap on the per-lattice loop.  In every
case the genus-level record is still written, so a guarded genus is bounded rather than
either running for hours or going missing.}
    genus_t0 := Realtime();
    data := Split(Split(Read(LabelPath("genera_basic", label)), "\n")[1], "|");
    basic_format := Split(Read("genera_basic.format"), "|");
    advanced_format := Split(Read("genera_advanced.format"), "|");
    hash_format := Split(Read("lat_hash.format"), "|");
    // This function only fills in basic lattice entries (essentially those that don't require interactions between different genera)
    lat_format := Split(Split(Read("lat_basic.format"), "\n")[1], "|");
    assert #data eq #basic_format;
    basics := AssociativeArray();
    for i in [1..#data] do
        basics[basic_format[i]] := data[i];
        if data[i] eq "None" then basics[basic_format[i]] := "\\N"; end if;
    end for;
    advanced := AssociativeArray();
    lats := [];

    n := StringToInteger(basics["rank"]);
    s := StringToInteger(basics["nplus"]);
    K := Rationals();
    LWG := LatticeWithGram;
    DualLat := Dual;
    rep := basics["rep"];
    // Switch to square brackets
    rep := "[" * rep[2..#rep - 1] * "]"; // Switch to square brackets
    gram0 := Matrix(K, n, eval rep);
    L0 := LWG(gram0 : CheckPositive := false);
    vprintf FillGenus, 1 : "Computing genus representatives...";
    reps := [];
    // Rank-2 lattices of square determinant -m^2 have isotropic (split) forms, on
    // which Magma's GenusRepresentatives fails (it reduces to binary quadratic
    // forms, which reject square discriminants); these are handled directly.  Note
    // the class number is NOT always 1 -- e.g. for m = 5 some genera have two
    // classes.
    // masslimit: a definite genus of mass M has class number on the order of M, so a
    // very large mass means enumeration is hopeless -- skip it and keep only the
    // (cheap) genus-level record.  Mass is only defined for definite genera; indefinite
    // ones store "\N" here, so only parse it when definite.
    skip_massive := false;
    if (masslimit gt 0) and (n eq s) then
        mass_pieces := Split(basics["mass"][2..#basics["mass"]-1], ",");
        genus_mass := StringToInteger(mass_pieces[1]) / StringToInteger(mass_pieces[2]);
        skip_massive := genus_mass gt masslimit;
    end if;
    if skip_massive then
        vprintf FillGenus, 1 : "Skipping enumeration: mass exceeds masslimit %o\n", masslimit;
        genus_success := false;
    elif n eq 2 and IsSquare(-Determinant(L0)) then
        genus_success, reps, elapsed := TimeoutCall(timeout, genus_reps_square_disc, <L0>, 1);
        vprintf FillGenus, 1 : "Genus representatives (square discriminant) computed in %o seconds\n", elapsed;
    elif (n eq s) and (n ge 3) then
        // Definite rank >= 3: our own p-neighbour enumeration (GenusRepresentativesFaster)
        // is the default.  It is mass-verified -- it enumerates until the accumulated
        // sum of 1/|Aut| equals the genus mass -- so it is provably complete when it
        // terminates.  It also avoids Magma's GenusRepresentatives, which is slow and
        // reliably fails at rank >= 7 (the "cs must be non-empty" bug); we used to run
        // it first and wait out its 60s timeout before falling back, ~60s wasted per
        // high-rank genus.  Magma is used only where Faster does not apply (below).
        genus_success, reps, elapsed := TimeoutCall(timeout, genus_reps_Faster, <L0>, 1);
        vprintf FillGenus, 1 : "Genus representatives (p-neighbours) computed in %o seconds\n", elapsed;
    elif (n ne s) and (n ge 3) then
        // Indefinite rank >= 3: SpinorRepresentatives is the default.  By Eichler each
        // proper spinor genus is a single class, so the spinor-genus representatives ARE
        // the genus representatives -- provably exact here.  It also avoids Magma's
        // GenusRepresentatives, which fails on some indefinite genera (e.g. signature
        // (2,2) rank 4 -- "Illegal null sequence" / "cs must be non-empty" in its LatNF
        // code) and is slower even when it works.
        genus_success, reps, elapsed := TimeoutCall(timeout, genus_reps_Spinor, <L0>, 1);
        vprintf FillGenus, 1 : "Genus representatives (spinor genera) computed in %o seconds\n", elapsed;
    else
        // Rank 1-2: neither Faster (needs definite rank >= 3) nor the spinor route (needs
        // rank >= 3, indefinite) applies, so use Magma's general GenusRepresentatives.
        genus_success, reps, elapsed := TimeoutCall(timeout, genus_reps_Magma, <L0>, 1);
        vprintf FillGenus, 1 : "Genus representatives computed in %o seconds\n", elapsed;
    end if;
    advanced["class_number"] := "\\N";
    advanced["adjacency_matrix"] := "\\N";
    advanced["adjacency_polynomials"] := "\\N";
    advanced["ambient_lattice"] := "\\N";   // TODO: compute the ambient lattice
    have_adjacency := false;   // set below iff the adjacency (Hecke) matrix is computed
    if genus_success then
        reps := reps[1];
        vprintf FillGenus, 1 : "Number of genus representatives: %o\n", #reps;
        advanced["class_number"] := #reps;
        vprintf FillGenus, 1 : "Computing adjacency matrix for p = ";
        hecke_mats := AssociativeArray();
        hecke_polys := AssociativeArray();
        G := Genus(L0);
        // This works for 2.28 - should be replaced by SetGenus in 2.29
        G`Representatives := reps;
        G`IsNatural := true;
        // The adjacency (Hecke) matrix is class_number x class_number and built by
        // p-neighbour walks -- an O(class_number^2) cost that dominates fill for
        // moderate class numbers at high rank (a single rank-8 class-~300 genus can take
        // hours).  Gate it on adjlimit, a much lower cap than sizelimit; the class number
        // is still recorded, and pneighbors below is emitted only when the matrix exists.
        have_adjacency := (n eq s) and ((adjlimit le 0) or (#reps le adjlimit));
        if have_adjacency then
          for p in hecke_primes(n) do
            vprintf FillGenus, 1 : "%o:", p;
            vtime FillGenus, 1 : Ap := AdjacencyMatrix(G,p);
            fpf := Factorization(CharacteristicPolynomial(Ap));
            hecke_mats[p] := Ap;
            hecke_polys[p] := [(<Coefficients(pair[1]), pair[2]>) : pair in fpf];
          end for;
          vprintf FillGenus, 1 : "Done!\n";
          advanced["adjacency_matrix"] := to_postgres(hecke_mats);
          advanced["adjacency_polynomials"] := to_postgres(hecke_polys);
        end if;
    else
        reps := [];
    end if;
    // sizelimit: past this class number, keep the genus-level record but skip the
    // per-lattice loop below (individual automorphism groups, theta series, hashes).
    if (sizelimit gt 0) and (#reps gt sizelimit) then
        vprintf FillGenus, 1 : "Class number %o exceeds sizelimit %o: storing genus-level data only\n", #reps, sizelimit;
        reps := [];
    end if;
    disc_invs := basics["discriminant_group_invs"];
    disc_invs := "[" * disc_invs[2..#disc_invs-1] * "]"; // Switch to square brackets
    disc_invs := eval disc_invs;
    disc_aut_size := #AutomorphismGroup(AbelianGroup(disc_invs));

    if (n eq s) then
        vprintf FillGenus, 1 : "Computing canonical forms and automorphism groups for representative ";
    end if;

    if (#reps gt 0) then
        to_per_rep := timeout div #reps + 1;
    end if;

    // Theta-series timing per precision, aggregated over the genus representatives.
    // DECISION: aggregate by max time per precision (the slowest representative
    // bounds the cost of hashing the whole genus to that precision).  This may
    // change in future -- e.g. to the mean or to per-representative timings.
    theta_elapsed := AssociativeArray();

    // timelimit: a per-genus wall-clock backstop.  If the per-lattice work overruns it
    // (a genus under sizelimit whose individual lattices are still pathologically slow),
    // stop and fall back to the genus-level record rather than grinding for hours.
    genus_timed_out := false;
    for Li->L in reps do
        if (timelimit gt 0) and (Realtime(genus_t0) gt timelimit) then
            vprintf FillGenus, 1 : "Per-genus timelimit %o s exceeded at rep %o/%o: storing genus-level data only\n", timelimit, Li, #reps;
            genus_timed_out := true;
            break;
        end if;
        lat := AssociativeArray();
        lat["lattice"] := L; // useful for subroutines; removed before saving to disk
        for col in ["rank", "nplus", "nminus", "disc_abs", "disc_sign", "disc_radical", "disc_witt", "disc_geometric", "disc_quadratic", "disc_half", "disc_2adic_unit", "bad_primes", "discriminant_group_invs", "discriminant_group_exponent", "is_even", "level", "scale", "conway_symbol", "dual_conway_symbol"] do
            lat[col] := basics[col];
        end for;
        lat["genus_label"] := basics["label"];
        lat["class_number"] := advanced["class_number"];
        // TODO (Eran) := The code for ConwaySymbol is currently in sage.
        // The magma implemntation is in version 2.29 that has some bugs
        // This is no longer part of the lattice, only of the genus
        // lat["dual_conway"] := "\\N";
        lat["aut_size"] := "\\N";
        lat["festi_veniani_index"] := "\\N";
        lat["aut_label"] := "\\N";
        lat["aut_group"] := "\\N";
        lat["is_chiral"] := "\\N";
        lat["orthogonal_complement"] := "\\N";
        lat["density"] := "\\N";
        lat["hermite"] := "\\N";
        lat["kissing"] := "\\N";
        lat["minimum"] := "\\N";
        lat["theta_series"] := "\\N";
        lat["theta_prec"] := "\\N";
        lat["dual_theta_series"] := "\\N";
        lat["dual_theta_prec"] := "\\N";
        lat["successive_minima"] := "\\N";
        // Only set below in the positive-definite (n eq s) branch; default them here
        // so the indefinite case still has every column.
        lat["center_density"] := "\\N";
        lat["canonical_gram"] := "\\N";
        lat["is_universal"] := "\\N";
        lat["is_even_universal"] := "\\N";
        // Trying to reduce the size of the entries in the gram matrix
        gram0 := GramMatrix(L);
        gram := LLLGram(gram0);
        max_abs := Max([Abs(x) : x in Eltseq(gram)]);
        max_abs_0 := Max([Abs(x) : x in Eltseq(gram0)]);
        if max_abs_0 le max_abs then
            gram := gram0;
        end if;
        lat["gram"] := Eltseq(gram);
        lat["gram_is_canonical"] := false;
        lat["gram_others"] := []; // This will be manually set in cases like E8 where we want to store other options
        // At the moment we do not have a notion of a canonical gram in the indefinite case
        if (n eq s) then
            vprintf FillGenus, 1 : "%o", gram;
            // The automorphism group and dual must be expressed in the same basis as
            // the stored gram (ConnectGenus rebuilds the lattice from that gram); once
            // the gram is canonicalised below, Lcanon holds the canonical lattice.
            Lcanon := L;
            success, canonical_gram, elapsed := TimeoutCall(to_per_rep, CanonicalForm, <gram>, 1);
            vprintf FillGenus, 1 : "Canonical form computed in %o seconds\n", elapsed;
            if success then
                canonical_gram := canonical_gram[1];
                lat["gram"] := Eltseq(canonical_gram);
                lat["canonical_gram"] := Eltseq(canonical_gram);
                lat["gram_is_canonical"] := true;
                Lcanon := LatticeWithGram(canonical_gram);
                lat["lattice"] := Lcanon;
            end if;
            success, aut_group, elapsed := TimeoutCall(to_per_rep, AutomorphismGroupFaster, <Lcanon>, 1);
            vprintf FillGenus, 1 : "Automorphism group computed in %o seconds\n", elapsed;
            if success then
                aut_group := aut_group[1];
                lat["aut_group"] := GroupToString(aut_group : use_id:=false);
                lat["aut_size"] := #aut_group;
                lat["is_chiral"] := &and[Determinant(g) eq 1 : g in Generators(aut_group)];
                // double checking, but also useful for festi-veniani
                LD := Dual(Lcanon : Rescale:=false);
                discL, quo := LD/Lcanon;
                disc_aut := AutomorphismGroup(discL);
                assert disc_aut_size eq #disc_aut;
                assert disc_invs eq Invariants(discL);
                gens_disc := [discL.i : i in [1..Ngens(discL)]];
                im_aut := [hom< discL -> discL | [quo(x@@quo*ChangeRing(alpha, Rationals())): x in gens_disc]> : alpha in Generators(aut_group)];
                lat["festi_veniani_index"] := disc_aut_size div #sub<disc_aut | im_aut>;
                if CanIdentifyGroup(#aut_group) then
                    Aid := IdentifyGroup(aut_group);
                    lat["aut_label"] := Sprintf("%o.%o", Aid[1], Aid[2]);
                end if;
            end if;
            lat["density"] := Density(L);
            lat["center_density"] := lat["density"] / SphereVolume(n);
            lat["hermite"] := HermiteNumber(L);
            lat["kissing"] := KissingNumber(L);
            m := Minimum(L);
            lat["minimum"] := m;
            target_prec := Max(150, m+4);
            theta, theta_prec, rep_theta_elapsed := ThetaSeriesIncremental(L, target_prec, to_per_rep);
            for p in Keys(rep_theta_elapsed) do
                theta_elapsed[p] := IsDefined(theta_elapsed, p)
                    // Maximum for worst case Hash computation.
                    // Change to sum if we care about the average case
                    select Maximum(theta_elapsed[p], rep_theta_elapsed[p]) else rep_theta_elapsed[p];
            end for;
            for p in Keys(theta_elapsed) do
                if not IsDefined(rep_theta_elapsed, p) then theta_elapsed[p] := Infinity(); end if;
            end for;
            lat["is_universal"] := "\\N";
            lat["is_even_universal"] := "\\N";
            if theta_prec gt 0 then
                lat["theta_series"] := theta;
                lat["theta_prec"] := theta_prec;
                // 15-theorem for is_universal, doubled 290-theorem for
                // is_even_universal (sources, and why 2*S290 is needed even for odd
                // lattices, in universality.m).  The previous scan certified
                // is_even_universal from theta support on 2..30, which is only valid
                // for doubly even lattices.  Norms beyond theta_prec are checked by
                // a TimeoutCall-guarded first-vector search; undecided stays "\N".
                udec, uval, eudec, euval := UniversalityFlags(L : Theta := theta);
                if udec then
                    lat["is_universal"] := uval;
                end if;
                if eudec then
                    lat["is_even_universal"] := euval;
                end if;
            else
                lat["theta_series"] := [1];
                lat["theta_prec"] := 1;
            end if;
            // Dual theta series (used as a tie-breaker in cmp_lat below, and stored
            // alongside the primal theta series).
            D := Dual(L);
            dual_theta, dual_theta_prec := ThetaSeriesIncremental(D, Max(150, Minimum(D) + 4), to_per_rep);
            if dual_theta_prec gt 0 then
                lat["dual_theta_series"] := dual_theta;
                lat["dual_theta_prec"] := dual_theta_prec;
            else
                lat["dual_theta_series"] := [1];
                lat["dual_theta_prec"] := 1;
            end if;
            //success, minima, elapsed := TimeoutCall(to_per_rep, SuccessiveMinima, <L>, 2);
            //vprintf FillGenus, 1 : "Successive minima computed in %o seconds\n", elapsed;
            //if success then 
            //lat["successive_minima"] := minima[1]; // For now, we throw away the vecs
            //end if;
            minima, vecs := SuccessiveMinima(L);
            lat["successive_minima"] := minima;

        end if;
        lat["hash"] := "\\N";

        //lat["level"] := Level(LatticeWithGram(ChangeRing(GramMatrix(L), Integers()) : CheckPositive:=false));

        Append(~lats, lat);
    end for;

    // If the per-genus timelimit tripped mid-loop, discard the partial lattices so the
    // record is consistent (class_number set, genus-level data only) rather than
    // claiming a class number it did not store.
    if genus_timed_out then lats := []; end if;

    vprintf FillGenus, 1 : "Done!\n";

    function cmp_lat(L1, L2)
        // 1. automorphism group size (larger first)
        if Type(L1["aut_size"]) eq RngIntElt and Type(L2["aut_size"]) eq RngIntElt then
            d := L2["aut_size"] - L1["aut_size"];
            if (d ne 0) then return d; end if;
        end if;
        // 2. density (denser first); within a genus the determinant is fixed, so we
        //    sort by the minimum as a proxy (larger minimum = denser)
        if Type(L1["minimum"]) eq RngIntElt and Type(L2["minimum"]) eq RngIntElt then
            d := L2["minimum"] - L1["minimum"];
            if (d ne 0) then return d; end if;
        end if;
        // 3. theta series
        if Type(L1["theta_series"]) eq SeqEnum and Type(L2["theta_series"]) eq SeqEnum then
            prec := Minimum(L1["theta_prec"], L2["theta_prec"]);
            for i in [1..prec - 1] do
                d := L1["theta_series"][i] - L2["theta_series"][i];
                if (d ne 0) then return d; end if;
            end for;
        end if;
        // 4. dual theta series
        if Type(L1["dual_theta_series"]) eq SeqEnum and Type(L2["dual_theta_series"]) eq SeqEnum then
            prec := Minimum(L1["dual_theta_prec"], L2["dual_theta_prec"]);
            for i in [1..prec - 1] do
                d := L1["dual_theta_series"][i] - L2["dual_theta_series"][i];
                if (d ne 0) then return d; end if;
            end for;
        end if;
        // 5. arbitrary tiebreaker: the (reduced) Gram matrix entries
        for i in [1..n^2] do
            d := L1["gram"][i] - L2["gram"][i];
            if (d ne 0) then return d; end if;
        end for;
        return 0;
    end function;

    // Tie breaker

    // Need dual_label, dual_conway
    // Compute festi_veniani_index in Sage?
    // Need label for lattice.  Don't want the label to rely on a difficult computation.  So we should probably avoid using the canonical form, and maybe avoid the automorphism group.
    // Proposal: Sort lexicographically by:
    // 1. Size of automorphism group (larger first): unfortunately this may be hard to compute
    // 2. Density - sort by minimum instead
    // 3. theta series
    // 4. dual theta series
    // 5. arbitrary tiebreaker
    // TODO: Sort reps according to canonical form?
    if (n eq s) then
        Sort(~lats, cmp_lat, ~perm);
    end if;

    SetColumns(0);
    for idx->L in lats do
        // Need label for lattice.
        lats[idx]["label"] := Sprintf("%o.%o", basics["label"], idx);
    end for;

    if #lats gt 0 then
        SetHashes(~lats, ~advanced, theta_elapsed, timeout);
    else
        // GenusRepresentatives failed (e.g. timed out on a hard genus): record the
        // genus with no stored lattices.  HashGenus factors through the canonical
        // label, so the genus hash is that of the label -- no representative needed.
        advanced["genus_hash"] := CollapseIntList(StringToBytes(basics["label"]));
        advanced["theta_distinguishing_prec"] := "\\N";
        advanced["is_theta_distinguished"] := "\\N";
        advanced["hash_function"] := "\\N";
        advanced["is_hash_distinguished"] := "\\N";
    end if;
    // We need to be able to look up hash functions for lattices that are not in the main
    // genus being processed.  So we write the hash function used to a separate file
    // so that it can be looked up when needed (see lookup_hash_function in connect_genus.m)
    Write(LabelPath("genera_hash", n, s, IntegerToString(advanced["genus_hash"]) : Create), advanced["hash_function"] : Overwrite);

    // TODO: Compute ambient_lattice

    for idx->L in lats do
        lat := L;
        if have_adjacency then   // pneighbors are read off the Hecke matrix, so only when it was computed
            pNeighbors := AssociativeArray();
            for p in hecke_primes(n) do
                // !!! TODO - check that permutation is applied in the right direction
                pNeighbors[p] := ["\"" * lats[j]["label"] * "\"" : j in [1..#lats] | hecke_mats[p][idx^perm,j^perm] ne 0];
            end for;
            lat["pneighbors"] := to_postgres(pNeighbors);
        else
            lat["pneighbors"] := "\\N";
        end if;
        Remove(~lat, "lattice");
        error if Keys(lat) ne Set(lat_format), [k : k in lat_format | k notin Keys(lat)], [k : k in Keys(lat) | k notin lat_format];
        output := Join([Sprintf("%o", to_postgres(lat[k])) : k in lat_format], "|");
        Write(LabelPath("lattice_basic_data", lat["label"] : Create), output : Overwrite);

    end for;
    // Now write hash data
    output := Join([Join([Sprintf("%o", to_postgres(lat[k])) : k in hash_format], "|") : lat in lats], "\n");
    Write(LabelPath("lattice_hashes", n, s, IntegerToString(advanced["genus_hash"]) : Create), output : Overwrite);

    error if Keys(basics) ne Set(basic_format), [k : k in basic_format | k notin Keys(basics)], [k : k in Keys(basics) | k notin basic_format];
    error if Keys(advanced) ne Set(advanced_format), [k : k in advanced_format | k notin Keys(advanced)], [k : k in Keys(advanced) | k notin advanced_format];
    output := Join([basics[k] : k in basic_format] cat [Sprintf("%o", advanced[k]) : k in advanced_format], "|");
    Write(LabelPath("genera_advanced", label : Create), output : Overwrite);
    return;
end intrinsic;
