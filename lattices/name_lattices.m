// Stage 4 of the pipeline: fill the "name" field of the lattices.
//
// Phase A names the "atomic" lattices that appear in Magma's LatticeDatabase
// (root lattices A_n/D_n/E_n and their duals, the laminated lattices LAMBDA_n,
// the Kappa lattices KAPPA_n / K12, Leech, Barnes-Wall, ...).  Phase B (to come)
// names the decomposable lattices from their orthogonal/tensor decompositions.

// A LatticeDatabase name is "clean" (a canonical name we want to use) if it is
// a family letter-block followed by a dimension and an optional dual star, e.g.
// A3, A3*, D4, E8, KAPPA8, LAMBDA12, BW16, Z2, I3, K12 -- or one of a few
// special names.  This deliberately rejects the database's variant spellings
// (e.g. "E8 (coding theory version)", "A6,1", "A5^{+2}", "LAMBDA11_MIN",
// "KAPPA14.1", "D4 as a Hurwitzian lattice").
function IsCleanName(s)
    if s in {"Leech", "Z"} then return true; end if;
    n := #s;
    if n eq 0 then return false; end if;
    if Substring(s, n, 1) eq "*" then n -:= 1; end if;   // drop a trailing dual star
    if n eq 0 then return false; end if;
    upper := "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    digits := "0123456789";
    i := 1;
    while i le n and Substring(s, i, 1) in upper do i +:= 1; end while;
    if i eq 1 or i gt n then return false; end if;   // need >=1 letter then >=1 digit
    while i le n do
        if not (Substring(s, i, 1) in digits) then return false; end if;
        i +:= 1;
    end while;
    return true;
end function;

intrinsic CleanNamedLattices() -> SeqEnum
{The canonically-named lattices in Magma's LatticeDatabase, as a sequence of
 tuples <name, lattice>, keeping only clean names (see IsCleanName).}//'
    D := LatticeDatabase();
    catalog := [];
    for i in [1..NumberOfLattices(D)] do
        nm := LatticeName(D, i);
        if IsCleanName(nm) then
            Append(~catalog, <nm, Lattice(D, i)>);
        end if;
    end for;
    return catalog;
end intrinsic;

// Family priority for resolving names that collide on the same lattice:
// Z > A > D > E > LAMBDA > KAPPA > special (smaller number = preferred).
function NameFamilyPriority(name)
    upper := "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    i := 1;
    while i le #name and Substring(name, i, 1) in upper do i +:= 1; end while;
    prefix := Substring(name, 1, i-1);
    case prefix:
        when "Z":      return 1;
        when "I":      return 2;
        when "A":      return 3;
        when "D":      return 4;
        when "E":      return 5;
        when "LAMBDA": return 6;
        when "KAPPA":  return 7;
        else           return 8;
    end case;
end function;

// The primitive integral form of L: clear denominators of the Gram matrix and
// divide out the content, giving the representative our database stores.
function PrimitiveIntegralForm(L)
    G := GramMatrix(L);
    d := LCM([Integers() | Denominator(x) : x in Eltseq(G)]);
    G := ChangeRing(d*G, Integers());
    G := G div GCD(Eltseq(G));
    return LatticeWithGram(G);
end function;

intrinsic AtomicName(L::Lat, catalog::SeqEnum) -> MonStgElt
{The canonical name of L if its primitive integral form is isometric to one of
 the named catalog lattices, resolving ties by family priority; "\N" otherwise.}
    Lp := PrimitiveIntegralForm(L);
    n := Rank(Lp);  det := Determinant(Lp);
    best := "\\N";  bestpri := 1000;
    for t in catalog do
        Cp := PrimitiveIntegralForm(t[2]);
        if Rank(Cp) ne n or Determinant(Cp) ne det then continue; end if;
        if IsIsometric(Lp, Cp) then
            pri := NameFamilyPriority(t[1]);
            if pri lt bestpri then best := t[1];  bestpri := pri; end if;
        end if;
    end for;
    return best;
end intrinsic;

intrinsic ScaledAtomicName(L::Lat, catalog::SeqEnum) -> MonStgElt
{The name of L when it is a scalar multiple of a named catalog lattice, with the
 scale written as a prefix, e.g. "A3", "A3*", "2A3"; "\N" if the primitive form
 of L is not a named lattice.}
    atomic := AtomicName(L, catalog);
    if atomic eq "\\N" then return "\\N"; end if;
    G := GramMatrix(L);
    if &and[ IsIntegral(x) : x in Eltseq(G) ] then
        c := GCD([Integers() | x : x in Eltseq(G)]);
        if c gt 1 then return Sprintf("%o%o", c, atomic); end if;
    end if;
    return atomic;
end intrinsic;

intrinsic NameAtomicLattices(DETCAP::RngIntElt : worker:=0, nworkers:=1) -> Assoc, Assoc
{The label -> name map for the named LatticeDatabase lattices -- and their integral
 scalings up to determinant DETCAP -- that occur in our database.  Each catalog
 lattice is located via its genus with FindLabel (much faster than testing every
 lattice against the catalog).  Name collisions are resolved by family priority;
 the second return value maps label -> priority (lower wins), which also lets the
 per-worker results be merged.  With nworkers > 1 the catalog is partitioned across
 workers (worker in [0..nworkers-1] does its share), so the naming can be split
 across parallel workers.}
    catalog := CleanNamedLattices();
    forms := [PrimitiveIntegralForm(t[2]) : t in catalog];
    // Balance the workers by the number of scalings each catalog entry generates
    // (~= (DETCAP/det)^(1/rank)); a rank-1 or low-det entry has vastly more than a
    // high-rank one, so a plain round-robin leaves some workers idle.  Every worker
    // sees the same catalog and runs the same deterministic greedy assignment.
    weights := [ Integers() | ];
    for i in [1..#catalog] do
        d := Determinant(forms[i]);  r := Rank(forms[i]);
        w := (d le DETCAP) select Floor(Root(RealField(30) ! DETCAP / d, r)) else 0;
        Append(~weights, w);
    end for;
    order := [ i : i in [1..#catalog] ];
    Sort(~order, func< a, b | weights[b] - weights[a] >);   // heaviest first
    wload := [ 0 : k in [1..nworkers] ];   // NB: `load` is a Magma keyword
    mine := {};
    for i in order do
        mw := 1;
        for k in [2..nworkers] do if wload[k] lt wload[mw] then mw := k; end if; end for;
        wload[mw] +:= weights[i];
        if mw eq worker + 1 then Include(~mine, i); end if;
    end for;

    names := AssociativeArray();
    prio  := AssociativeArray();
    for ci in mine do
        name := catalog[ci][1];
        Lp := forms[ci];
        r := Rank(Lp);  d := Determinant(Lp);  p := NameFamilyPriority(name);
        c := 1;
        while d * c^r le DETCAP do
            L := c eq 1 select Lp else LatticeWithGram(c * GramMatrix(Lp));
            label := FindLabel(L);
            newname := c gt 1 select Sprintf("%o%o", c, name) else name;
            // Resolve collisions by (priority, name): lower priority wins, ties go to
            // the lexicographically smaller name (e.g. "A2" over its dual "A2*").  This
            // is order-independent, so per-worker results merge to the serial answer.
            if label ne "\\N" and ((not IsDefined(prio, label)) or p lt prio[label]
                                   or (p eq prio[label] and newname lt names[label])) then
                names[label] := newname;
                prio[label] := p;
            end if;
            c +:= 1;
        end while;
    end for;
    return names, prio;
end intrinsic;

// --- Phase B: composing names of decomposable lattices ---------------------
// Operator precedence (tightest first): scale prefix / *, then ^, then x, then +.

// Parenthesise a name before raising it to a power (^): needed if it already
// contains a lower-or-equal precedence operator (+ or x).
function WrapForPower(s)
    return ("+" in s) or ("x" in s) select "(" cat s cat ")" else s;
end function;

// Parenthesise a name as a tensor factor (x): needed only if it contains a sum.
function WrapForTensor(s)
    return ("+" in s) select "(" cat s cat ")" else s;
end function;

intrinsic ComposeSumName(factors::SeqEnum) -> MonStgElt
{The orthogonal-sum name from ordered factors given as tuples <name, multiplicity>
 (caller orders by descending dimension, then family priority).  Repeated factors
 use power notation, e.g. <"D4",1>,<"A2",2> -> "D4+A2^2".}
    terms := [ f[2] gt 1 select Sprintf("%o^%o", WrapForPower(f[1]), f[2]) else f[1]
               : f in factors ];
    return Join(terms, "+");
end intrinsic;

intrinsic ComposeTensorName(factors::SeqEnum) -> MonStgElt
{The tensor-product name from ordered factor names, parenthesising any summand
 factors, e.g. ["A2","A2"] -> "A2xA2", ["A2+D4","A2"] -> "(A2+D4)xA2".}
    return Join([ WrapForTensor(nm) : nm in factors ], "x");
end intrinsic;

// Order factor labels by descending dimension (the label's rank field), then by
// family priority of the assigned name.  Returns the permutation of indices.
function FactorOrder(labels, names)
    order := [1..#labels];
    Sort(~order, func< i, j |
        di ne dj select dj - di else NameFamilyPriority(names[labels[i]]) - NameFamilyPriority(names[labels[j]])
            where di := StringToInteger(Split(labels[i], ".")[1])
            where dj := StringToInteger(Split(labels[j], ".")[1]) >);
    return order;
end function;

intrinsic OrthogonalSumName(factor_labels::SeqEnum, mults::SeqEnum, names::Assoc) -> MonStgElt
{Compose the orthogonal-sum name of a lattice from the labels of its indecomposable
 factors and their multiplicities, using the label->name map.  Factors are ordered
 by descending dimension then family priority.  Returns "\N" if some factor is not
 yet named.}
    if exists{ f : f in factor_labels | not IsDefined(names, f) } then return "\\N"; end if;
    order := FactorOrder(factor_labels, names);
    return ComposeSumName([ <names[factor_labels[i]], mults[i]> : i in order ]);
end intrinsic;

// Split s on a separator character, ignoring separators nested inside [ ].
function SplitDepth0(s, sep)
    parts := [];  depth := 0;  start := 1;
    for i in [1..#s] do
        c := Substring(s, i, 1);
        if c eq "[" then depth +:= 1;
        elif c eq "]" then depth -:= 1;
        elif c eq sep and depth eq 0 then
            Append(~parts, Substring(s, start, i-start));  start := i+1;
        end if;
    end for;
    Append(~parts, Substring(s, start, #s-start+1));
    return parts;
end function;

intrinsic ParseTensorDecompositions(s::MonStgElt) -> SeqEnum
{Parse the tensor_decompositions field "[[["lab",ct],...],...]" (labels unquoted
 before 2026-11) into a sequence of
 options, each a sequence of tuples <factor_label, count>; "\N" -> empty.}
    if s eq "\\N" or #s lt 2 then return []; end if;
    result := [];
    for optstr in SplitDepth0(Substring(s, 2, #s-2), ",") do
        opt := [];
        for pairstr in SplitDepth0(Substring(optstr, 2, #optstr-2), ",") do
            lc := Split(Substring(pairstr, 2, #pairstr-2), ",");
            lab := lc[1];
            if #lab ge 2 and lab[1] eq "\"" then lab := Substring(lab, 2, #lab-2); end if;   // quoted since 2026-11 (JSON)
            Append(~opt, <lab, StringToInteger(lc[2])>);
        end for;
        Append(~result, opt);
    end for;
    return result;
end intrinsic;

intrinsic TensorProductName(option::SeqEnum, names::Assoc) -> MonStgElt
{Compose the tensor-product name from one tensor-decomposition option, given as a
 sequence of tuples <factor_label, count>, using the label->name map.  Factors are
 ordered by descending dimension then family priority and expanded by count.
 Returns "\N" if some factor is not yet named.}
    labels := [ p[1] : p in option ];
    if exists{ f : f in labels | not IsDefined(names, f) } then return "\\N"; end if;
    counts := AssociativeArray();
    for p in option do counts[p[1]] := p[2]; end for;
    order := FactorOrder(labels, names);
    factor_names := [];
    for i in order do
        for k in [1 .. counts[labels[i]]] do Append(~factor_names, names[labels[i]]); end for;
    end for;
    return ComposeTensorName(factor_names);
end intrinsic;

// --- naming a single lattice during ConnectGenus (stage 4 composite part) ---

intrinsic LoadAtomicNames() -> Assoc
{Load the atomic-name map written by run_basic_names (the file "atomic_names",
 one "label|name" per line) as a label -> name associative array.}
    names := AssociativeArray();
    if OpenTest("atomic_names", "r") then
        for line in Split(Read("atomic_names"), "\n") do
            if #line eq 0 then continue; end if;
            parts := Split(line, "|");
            names[parts[1]] := parts[2];
        end for;
    end if;
    return names;
end intrinsic;

intrinsic LatticeName(label::MonStgElt, ortho_factors::SeqEnum, ortho_mults::SeqEnum,
                      atomic_names::Assoc, name_i::RngIntElt) -> MonStgElt
{The name of the lattice with the given label: its atomic name if present in the
 atomic_names map, otherwise composed from its decompositions (tensor product
 preferred over orthogonal sum), reading each factor's name from atomic_names or,
 failing that, from the factor's lattice_advanced_data file (column name_i).
 Returns "\N" if it cannot yet be named (some factor is not yet named).}
    if IsDefined(atomic_names, label) then return atomic_names[label]; end if;

    factor_name := function(f)
        if f eq "\\N" then return "\\N"; end if;   // a factor not in the database
        if IsDefined(atomic_names, f) then return atomic_names[f]; end if;
        path := LabelPath("lattice_advanced_data", f);
        if not OpenTest(path, "r") then return "\\N"; end if;
        return Split(Split(Read(path), "\n")[1], "|")[name_i];
    end function;

    // tensor products preferred
    pieces := Split(label, ".");
    tpath := "lattice_decomp_data/" * pieces[1] * "." * pieces[2];
    if OpenTest(tpath, "r") then
        for line in Split(Read(tpath), "\n") do
            if #line eq 0 then continue; end if;
            parts := Split(line, "|");
            if parts[1] ne label then continue; end if;
            for opt in ParseTensorDecompositions(parts[2]) do
                fnames := AssociativeArray();
                for pr in opt do
                    nm := factor_name(pr[1]);
                    if nm ne "\\N" then fnames[pr[1]] := nm; end if;
                end for;
                t := TensorProductName(opt, fnames);
                if t ne "\\N" then return t; end if;
            end for;
            break;
        end for;
    end if;

    // orthogonal sum
    if #ortho_factors gt 0 then
        fnames := AssociativeArray();
        for f in ortho_factors do
            nm := factor_name(f);
            if nm ne "\\N" then fnames[f] := nm; end if;
        end for;
        return OrthogonalSumName(ortho_factors, ortho_mults, fnames);
    end if;

    return "\\N";
end intrinsic;
