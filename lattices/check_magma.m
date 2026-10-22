// Refuse to run the census on a Magma whose GenusRepresentatives is wrong on
// indefinite binary lattices (see MagmaBinaryGenusPatched in genus_reps.m).
// Usage: magma -b check_magma.m.  Patched iff it prints MAGMA_OK: an error can
// still exit 0 under SetQuitOnError, so callers must not trust the exit status alone.
SetQuitOnError(true);
AttachSpec("lattices.spec");
a, b, c := GetVersion();
version := Sprintf("%o.%o-%o", a, b, c);
ok, msg := MagmaBinaryGenusPatched();
if not ok then
    printf "ERROR: this Magma (%o) is not patched: %o\nApply magma_patches/neighbors-2.30-1.patch to its package directory.\n", version, msg;
    exit 1;
end if;
printf "MAGMA_OK: Magma %o gets indefinite binary genera right\n", version;
exit 0;
