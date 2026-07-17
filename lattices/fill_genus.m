declare verbose FillGenus, 1;

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

function adjacency_matrix_of(G, p)
    // Wrapped so it can run under TimeoutCall: AdjacencyMatrix walks the p-neighbour graph
    // in a C-level builtin that does not yield to Magma's Alarm, so it needs the HardKill
    // watchdog to be bounded (it can take >75s for a single rank-8 genus).
    return AdjacencyMatrix(G, p);
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

intrinsic FillGenus(label::MonStgElt : timeout := 1800, masslimit := 0, sizelimit := 0, timelimit := 0, adjlimit := 0, perlattice := 0, use_orth := true, force_orth := false)
{Fill the data for a genus and its lattice representatives, given files in the genera_basic format.

Enumeration guards (0 = unlimited): masslimit skips enumerating a definite genus whose
mass exceeds it (its class number is too large to enumerate); adjlimit is a work budget
for the adjacency (Hecke) matrix -- skip it when the estimated work, ~class_number times
sum_p p^(rank-2) over the Hecke primes p, exceeds adjlimit;
sizelimit records the genus-level data but skips storing individual lattices past that
class number; timelimit is a per-genus wall-clock cap on the per-lattice loop.  In every
case the genus-level record is still written, so a guarded genus is bounded rather than
either running for hours or going missing.

perlattice is the wall-clock each individual lattice gets for its own computations
(canonical form, automorphism group, theta series).  It used to be derived as
timeout/class_number, which starved exactly the genera that needed it most -- at
timeout 300 a class-number-55 genus gave each lattice 6s, while CanonicalForm alone
measures 8-66s at rank 12, so those fields were silently written as "\N".  Note this
trades against timelimit: the per-lattice loop costs roughly class_number*perlattice,
and a genus that overruns timelimit has ALL of its lattices discarded, so raising
perlattice without raising timelimit turns partial data into no data.  Left at 0 it
falls back to the old timeout/class_number budget.

use_orth selects the orbit-method
enumeration for definite genera of rank >= 3 (see GenusReps in genus_reps.m).}
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
    // masslimit: a definite genus of mass M has class number on the order of M, so a
    // very large mass means enumeration is hopeless -- skip it and keep only the
    // (cheap) genus-level record.  Mass is only defined for definite genera; indefinite
    // ones store "\N" here, so only parse it when definite.
    genus_mass := 0;
    if (n eq s) and (basics["mass"] ne "\\N") then
        mass_pieces := Split(basics["mass"][2..#basics["mass"]-1], ",");
        genus_mass := StringToInteger(mass_pieces[1]) / StringToInteger(mass_pieces[2]);
    end if;
    skip_massive := (masslimit gt 0) and (n eq s) and (genus_mass gt masslimit);
    if skip_massive then
        vprintf FillGenus, 1 : "Skipping enumeration: mass exceeds masslimit %o\n", masslimit;
        genus_success := false;
    else
        // Method selection, per-method timing and fallbacks live in GenusReps
        // (genus_reps.m); use_orth toggles the orbit method for definite rank >= 3,
        // refined by the rank/mass heuristic (the neighbour walk stays the default
        // for small class numbers at low rank, where it wins).
        genus_success, reps, orth_prov := GenusReps(L0 : Timeout := timeout,
            UseOrth := use_orth and (force_orth or UseOrthHeuristic(n, genus_mass)));
    end if;
    advanced["class_number"] := "\\N";
    advanced["adjacency_matrix"] := "\\N";
    advanced["adjacency_polynomials"] := "\\N";
    // ambient_genus: the genus of the parent lattices from which this genus's
    // classes were obtained as orthogonal complements (orbit-method provenance);
    // the per-lattice counter/vector references live in lat_lattices.
    if not assigned orth_prov then orth_prov := 0; end if;
    advanced["ambient_genus"] := (orth_prov cmpeq 0 or orth_prov[1] eq "") select "\\N" else orth_prov[1];
    have_adjacency := false;   // set below iff the adjacency (Hecke) matrix is computed
    if genus_success then
        vprintf FillGenus, 1 : "Number of genus representatives: %o\n", #reps;
        advanced["class_number"] := #reps;
        vprintf FillGenus, 1 : "Computing adjacency matrix for p = ";
        hecke_mats := AssociativeArray();
        hecke_polys := AssociativeArray();
        G := Genus(L0);
        // This works for 2.28 - should be replaced by SetGenus in 2.29
        G`Representatives := reps;
        G`IsNatural := true;
        // The adjacency (Hecke) matrix walks the p-neighbour graph: each of the #reps
        // representatives has ~p^(n-2) p-neighbours (the isotropic lines in L/pL), each
        // needing an isometry test.  So the work is ~ #reps * sum_p p^(n-2), which grows
        // as p^(n-2) with rank -- a single rank-12 genus costs 16x a rank-8 one even at
        // class number 1.  We estimate that work and skip the matrix past adjlimit (a
        // work budget), an adaptive cutoff that captures both the rank and class-number
        // dependence, rather than a crude cap on either alone.  The class number is still
        // recorded, and pneighbors below is emitted only when the matrix exists.
        adj_work := #reps * (&+[Integers() | p^Maximum(n-2, 0) : p in hecke_primes(n)]);
        have_adjacency := (n eq s) and ((adjlimit le 0) or (adj_work le adjlimit));
        if have_adjacency then
          adj_ok := true;
          for p in hecke_primes(n) do
            vprintf FillGenus, 1 : "%o:", p;
            // Bound the adjacency computation: it is an uninterruptible builtin that can
            // run for well over a minute on a single high-rank genus, so it needs the
            // HardKill watchdog.  If it overruns we drop the (optional) Hecke data for
            // this genus -- the class number is already recorded -- rather than let a few
            // genera dominate wall-clock.
            adj_success, adj_out, adj_elapsed := TimeoutCall(timeout, adjacency_matrix_of, <G, p>, 1 : HardKill := true);
            if not adj_success then
                vprintf FillGenus, 1 : "(timed out after ~%os, skipping Hecke data)\n", timeout;
                adj_ok := false;
                break;
            end if;
            Ap := adj_out[1];
            vprintf FillGenus, 1 : "computed in %o seconds\n", adj_elapsed;
            fpf := Factorization(CharacteristicPolynomial(Ap));
            hecke_mats[p] := Ap;
            hecke_polys[p] := [(<Coefficients(pair[1]), pair[2]>) : pair in fpf];
          end for;
          if adj_ok then
            vprintf FillGenus, 1 : "Done!\n";
            advanced["adjacency_matrix"] := to_postgres(hecke_mats);
            advanced["adjacency_polynomials"] := to_postgres(hecke_polys);
          else
            have_adjacency := false;   // no Hecke data -> the pneighbors block below is skipped
          end if;
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
        // PER-OPERATION timeout, NOT a per-lattice aggregate: to_per_rep is passed
        // independently to CanonicalForm, AutomorphismGroupFaster, the primal theta series
        // and the dual theta series below, so a single lattice can spend up to ~4*to_per_rep
        // (plus the unbounded Minimum/kissing work between them).  The aggregate wall-clock
        // guard is the per-GENUS timelimit checked at the top of this loop, not this value.
        // perlattice gives every operation the same budget regardless of class number; the
        // old timeout/#reps split a fixed budget across the representatives, which starved
        // precisely the large genera (class number 55 at timeout 300 => 6s per operation,
        // against a measured 8-66s just for CanonicalForm at rank 12) and silently wrote the
        // results as "\N".  Kept as the fallback when perlattice is 0.
        to_per_rep := (perlattice gt 0) select perlattice else timeout div #reps + 1;
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
        // Orbit-method provenance: this lattice is the orthogonal complement of
        // orthogonal_complement (a vector) inside the lattice with counter
        // ambient_lattice in the genus ambient_genus (see genera_advanced).  The
        // vector lives in the first lattice up the parent chain for which a Gram
        // matrix is available; \N when the class came from a p-neighbour walk or
        // the parent genus has no recorded label/counters.
        lat["ambient_lattice"] := "\\N";
        lat["orthogonal_complement"] := "\\N";
        if orth_prov cmpne 0 and orth_prov[1] ne "" and Li le #orth_prov[2] and orth_prov[2][Li][1] gt 0 then
            lat["ambient_lattice"] := orth_prov[2][Li][1];
            lat["orthogonal_complement"] := orth_prov[2][Li][2];
        end if;
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
                if lat["is_even"] eq "T" then
                    lat["is_universal"] := false;
                    lat["is_even_universal"] := "\\N";
                    // 15 theorem
                    for m in [2..Min(30,#theta-1) by 2] do
                        if theta[m+1] eq 0 then
                            lat["is_even_universal"] := false;
                            break;
                        end if;
                    end for;
                    if lat["is_even_universal"] cmpeq "\\N" and #theta ge 31 then
                        lat["is_even_universal"] := true;
                    end if;
                else
                    // 290 theorem
                    for m in [1..Min(290,#theta-1)] do
                        if theta[m+1] eq 0 then
                            lat["is_universal"] := false;
                            break;
                        end if;
                    end for;
                    if lat["is_universal"] cmpeq "\\N" and #theta ge 291 then
                        lat["is_universal"] := true;
                    end if;
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
        // Need label for lattice.  counter records this sort position explicitly, so
        // cross-genus references (ambient_lattice in child genera) can point at it.
        lats[idx]["label"] := Sprintf("%o.%o", basics["label"], idx);
        lats[idx]["counter"] := idx;
    end for;
    // Persist the finished genus (reps in counter order) to the on-disk orth cache:
    // in rank-descending orchestration the genera of the next-lower rank descend
    // from these lattices and resolve their ambient_lattice counters against them.
    if (n eq s) and (#lats gt 0) then
        WriteGenusOrthCache(basics["label"], [lat["lattice"] : lat in lats]);
    end if;

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
