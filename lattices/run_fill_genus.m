// This script is used to find all of the representatives in a positive definite genus, along with some basic quantities of each lattice in the genus (roughly those that do not require access to lattices from any other genus) and some additional attributes of the genus itself.
// Usage: magma -b label:=foo run_fill_genus.m
// Batch: magma -b labels:=foo:bar:baz run_fill_genus.m
// Options: timeout:=N (default 60, seconds per label)
//
// Parallel across servers:
//   xargs -n 100 < genera_todo.txt | tr ' ' ':' > genera_todo_chunked.txt
//   parallel --sshloginfile servers.txt --joblog jobs.log --eta --resume \
//     'cd ~/projects/k3s-lmfdb/lattices && magma -b labels:={} verbose:=0 run_fill_genus.m' \
//     :::: genera_todo_chunked.txt > output.txt
//
// Errors are prefixed with "ERROR: label: ..."
// Extract retry list: grep '^ERROR:' output.txt | cut -d: -f2 | tr -d ' ' > genera_failed.txt

AttachSpec("lattices.spec");

if not assigned verbose then verbose := "0"; end if;
SetVerbose("FillGenus", StringToInteger(verbose));
SetVerbose("GenusReps", StringToInteger(verbose));

if not assigned timeout then timeout := "60"; end if;
timeout := StringToInteger(timeout);

// Enumeration guards (0 = unlimited): see FillGenus.
if not assigned masslimit then masslimit := "0"; end if;
masslimit := StringToInteger(masslimit);
if not assigned sizelimit then sizelimit := "0"; end if;
sizelimit := StringToInteger(sizelimit);
if not assigned timelimit then timelimit := "0"; end if;
timelimit := StringToInteger(timelimit);
if not assigned adjlimit then adjlimit := "0"; end if;
adjlimit := StringToInteger(adjlimit);
// Wall-clock each individual lattice gets (canonical form, automorphism group, theta).
// 0 falls back to the old timeout/class_number budget.  See FillGenus.
if not assigned perlattice then perlattice := "0"; end if;
perlattice := StringToInteger(perlattice);
// useorth:=0 disables the orbit-method enumeration; 1 (default) routes by
// UseOrthHeuristic; 2 forces the orbit method for every definite rank >= 3 genus
// (mainly for testing/benchmarking, and to exercise the provenance columns at
// ranks the heuristic would route to the neighbour walk).
if not assigned useorth then useorth := "1"; end if;
useorth := StringToInteger(useorth);
forceorth := useorth ge 2;
useorth := useorth ge 1;
// warmlimit: wall-clock budget (seconds) for the batch cache pre-warm below.
if not assigned warmlimit then warmlimit := "300"; end if;
warmlimit := StringToInteger(warmlimit);

if assigned labels then
    label_list := Split(labels, ":");
else
    label_list := [label];
end if;

procedure() // forcing magma to read the full input before forking
// Pre-warm the orth cache across the whole batch before the per-label loop:
// one descent sweep can serve several genera at once, and each FillGenus call
// runs its enumeration inside a TimeoutCall fork, which inherits this
// process's warm cache.  Failures here are non-fatal -- FillGenus enumerates
// on its own if the cache is cold.
if useorth and #label_list gt 1 then
    // The pre-warm cannot run under TimeoutCall (the fork would take the warmed
    // cache with it) and Alarm is uncatchable, so it is bounded by an internal
    // wall-clock budget instead: past warmlimit seconds the orth machinery raises
    // a catchable error at its next safe point and we proceed with whatever part
    // of the cache was already warmed (plus the disk tier).
    SetGenusOrthTimeBudget(warmlimit);
    try
        basic_format := Split(Read("genera_basic.format"), "|");
        ranki := Index(basic_format, "rank");
        nplusi := Index(basic_format, "nplus");
        repi := Index(basic_format, "rep");
        massi := Index(basic_format, "mass");
        batch := [];
        for l in label_list do
            data := Split(Split(Read(LabelPath("genera_basic", l)), "\n")[1], "|");
            n := StringToInteger(data[ranki]);
            if data[nplusi] ne data[ranki] or n lt 3 then continue; end if;  // definite rank >= 3 only
            gm := 0;
            mstr := data[massi];
            if mstr ne "\\N" and mstr ne "None" and #mstr gt 2 then
                mp := Split(mstr[2..#mstr-1], ",");
                gm := StringToInteger(mp[1]) / StringToInteger(mp[2]);
            end if;
            // only pre-warm genera that the dispatch will route to the orbit method
            if not (forceorth or UseOrthHeuristic(n, gm)) then continue; end if;
            Append(~batch, GramStringToLat(data[repi], n));
        end for;
        if #batch gt 1 then
            _ := GenusRepresentativesOrthBatch(batch);
        end if;
    catch e
        printf "WARNING: batch pre-warm stopped early (%o); continuing with a partial cache\n", e`Object;
    end try;
    ClearGenusOrthTimeBudget();
end if;
failed := [];
for l in label_list do
    try
        FillGenus(l : timeout := timeout, masslimit := masslimit, sizelimit := sizelimit, timelimit := timelimit, adjlimit := adjlimit, perlattice := perlattice, use_orth := useorth, force_orth := forceorth);
    catch e
        printf "ERROR: %o: %o\n", l, e;
        Append(~failed, l);
    end try;
end for;

if #failed gt 0 then
    exit 1;
end if;
exit 0;
end procedure();
