// Refuse to run the census on a Magma whose GenusRepresentatives is wrong on
// indefinite binary lattices, whose Mass is wrong in ranks 2 mod 4, or whose spinor
// characters are wrong on indefinite lattices of rank >= 3 (MagmaBinaryGenusPatched,
// MagmaMassPatched, MagmaSpinorPatched in genus_reps.m).
// Usage: magma -b check_magma.m.  Patched iff it prints MAGMA_OK: an error can
// still exit 0 under SetQuitOnError, so callers must not trust the exit status alone.
SetQuitOnError(true);
AttachSpec("lattices.spec");
procedure()  // parse the whole file before GenusReps forks (CLAUDE.md, TimeoutCall gotcha)
a, b, c := GetVersion();
version := Sprintf("%o.%o-%o", a, b, c);
// First, before this process loads the patched packages itself: enumerate a definite
// genus exactly as fill does (GenusReps forks a HardKill child).  A patched package whose
// .sig is older than its .m can make Magma recompile it on first use, which needs write
// access to the Magma install.  On lovelace (2026-11-04) Mass, SpinorCharacters and
// GenusRepresentatives loaded fine, in this process and in a child, but fill's
// GenusReps could not load neighbors.m, so every Magma-enumerated definite genus of
// the calibration was enum_fail.
ok, reps := GenusReps(StandardLattice(9) : Timeout := 600, UseOrth := false);
if not ok or #reps ne 2 then
    printf "ERROR: this Magma (%o) cannot enumerate I_9 the way fill does (GenusReps: %o, %o classes, expected 2).  If the output above says it cannot write a .lck file, the owner of the Magma install must recompile the patched packages (see CLAUDE.md).\n", version, ok, #reps;
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
end procedure();
quit;
