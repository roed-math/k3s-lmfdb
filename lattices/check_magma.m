// Refuse to run the census on a Magma whose GenusRepresentatives is wrong on
// indefinite binary lattices, whose Mass is wrong in ranks 2 mod 4, or whose spinor
// characters are wrong on indefinite lattices of rank >= 3 (MagmaBinaryGenusPatched,
// MagmaMassPatched, MagmaSpinorPatched in genus_reps.m).
// Usage: magma -b check_magma.m.  Patched iff it prints MAGMA_OK: an error can
// still exit 0 under SetQuitOnError, so callers must not trust the exit status alone.
SetQuitOnError(true);
AttachSpec("lattices.spec");
a, b, c := GetVersion();
version := Sprintf("%o.%o-%o", a, b, c);
// First, before this process loads the patched packages itself: use them in a HardKill
// child, as fill does (GenusReps forks before enumerating).  A patched package whose
// .sig is older than its .m makes Magma recompile it on first use, which needs write
// access to the Magma install; on lovelace (2026-11-04) the parent loaded neighbors.m
// fine, but fill's child could not, so every Magma-enumerated definite genus failed.
probe := function(L)
    _ := Mass(StandardLattice(14));
    T := LatticeWithGram(Matrix(Integers(), 3, 3, [4,1,0, 1,10,0, 0,0,-54]) : CheckPositive := false);
    _ := SpinorCharacters(Genus(T) : Proper := false);
    return #GenusRepresentatives(L);
end function;
ok, out := TimeoutCall(600, probe, <StandardLattice(9)>, 1 : HardKill := true);
if not ok or out[1] ne 2 then
    printf "ERROR: this Magma (%o) cannot use its lattice packages in a child process (got %o): if the output above says it cannot write a .lck file, have the owner of the Magma install run this script once so the patched packages are recompiled.\n", version, ok select out[1] else "a failure";
    exit 1;
end if;
ok, msg := MagmaBinaryGenusPatched();
if not ok then
    printf "ERROR: this Magma (%o) is not patched: %o\nApply magma_patches/neighbors-2.30-1.patch to its package directory.\n", version, msg;
    exit 1;
end if;
ok, msg := MagmaMassPatched();
if not ok then
    printf "ERROR: this Magma (%o) has the Mass bug: %o\nApply magma_patches/mass-2.30-1.patch to its package directory.\n", version, msg;
    exit 1;
end if;
ok, msg := MagmaSpinorPatched();
if not ok then
    printf "ERROR: this Magma (%o) has the spinor character bugs: %o\nApply magma_patches/diagonal-2.30-1.patch to its package directory.\n", version, msg;
    exit 1;
end if;
printf "MAGMA_OK: Magma %o has the binary genus, mass and spinor fixes\n", version;
exit 0;
