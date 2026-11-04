// The theta-Hecke pass: decide the genera whose theta series fill left undecided
// (SetHashes writes theta_hecke/<rank>/<nplus>/<label> for them), with the p-neighbour
// graph and the Hecke eigenspaces of the cusp space (theta_hecke.m).  Run from lattices/
// after fill; run_all.py groups the genera by cusp space, one file per job.
//
// Usage: magma -b labels:=FILE [timeout:=600] [memlimit:=32] run_theta_hecke.m
//   FILE: genus labels, one per line; genera sharing a space should be in one file
//   timeout: seconds per adjacency computation and per cusp space (HardKill)
//   memlimit: GB for this process and its children, 0 = none
// Writes theta_hecke_done/<rank>/<nplus>/<label> (see ThetaHeckePass); errors are
// printed as "ERROR: <file>: ..." and exit 1.

AttachSpec("lattices.spec");
if not assigned memlimit then memlimit := "32"; end if;
if StringToInteger(memlimit) gt 0 then SetMemoryLimit(StringToInteger(memlimit) * 2^30); end if;
if not assigned timeout then timeout := "600"; end if;

procedure()   // read the whole script before the first fork
    label_list := [l : l in Split(Read(labels), "\n") | l ne ""];
    try
        ThetaHeckePass(label_list : timeout := StringToInteger(timeout));
    catch e
        printf "ERROR: %o: %o\n", labels, e;
        quit 1;
    end try;
end procedure();
quit;
