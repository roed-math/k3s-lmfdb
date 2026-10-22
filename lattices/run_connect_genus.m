// This script is run after run_fill_genus.m, and fills in additional attributes of lattices in the genus that need access to lattices in different genera.
// Usage: magma -b label:=foo run_connect_genus.m
// Batch: magma -b labels:=foo:bar:baz run_connect_genus.m
// Options: timeout:=N (default 60, seconds per label)
//
// Parallel across servers:
//   xargs -n 100 < genera_todo.txt | tr ' ' ':' > genera_todo_chunked.txt
//   parallel --sshloginfile servers.txt --joblog jobs.log --eta --resume \
//     'cd ~/projects/k3s-lmfdb/lattices && magma -b labels:={} verbose:=0 run_connect_genus.m' \
//     :::: genera_todo_chunked.txt > output.txt
//
// Errors are prefixed with "ERROR: label: ..."
// Extract retry list: grep '^ERROR:' output.txt | cut -d: -f2 | tr -d ' ' > genera_failed.txt

AttachSpec("lattices.spec");

// memlimit:=G caps this process (and every TimeoutCall child it forks) at G GB of
// memory; 0 = no cap.  A child over the cap dies and its computation counts as
// failed; the parent over the cap ends the batch, which run_fill_genus.sh reports.
if not assigned memlimit then memlimit := "0"; end if;
if StringToInteger(memlimit) gt 0 then SetMemoryLimit(StringToInteger(memlimit) * 2^30); end if;
if not assigned verbose then verbose := "0"; end if;
SetVerbose("ConnectGenus", StringToInteger(verbose));

if not assigned timeout then timeout := "60"; end if;
timeout := StringToInteger(timeout);

// phase controls the two-pass connect (see run_all.py): "1" = produce (compute
// everything except the decomposable-lattice derivations, writing the indecomposable
// factor data); "2" = consume (derive the decomposable fields from factor data now on
// disk); unset/"0" = single pass (do everything).
if not assigned phase then phase := "0"; end if;

if assigned labels then
    label_list := Split(labels, ":");
else
    label_list := [label];
end if;

procedure() // forcing magma to read the full input before forking
failed := [];
for l in label_list do
    try
        if phase eq "2" then
            ConnectGenusDerive(l);
        elif phase eq "1" then
            ConnectGenus(l : timeout := timeout, derive_decomposable := false);
        else
            ConnectGenus(l : timeout := timeout);
        end if;
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
