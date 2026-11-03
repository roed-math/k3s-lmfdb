# Fixes to Hecke applied at load time (from __init__, in the master and in every pool
# worker), until a Hecke release that contains them is pinned in oscar/Manifest.toml.
#
# Hecke.jl PR #2611 (https://github.com/thofma/Hecke.jl/pull/2611, open on 2026-11-03):
# `_automorphous_numbers` added the image of -1 in each local factor instead of its
# diagonal image, and `automorphous_numbers` sorted a copy of I instead of I.  With the
# pinned Hecke 0.39.21, `improper_spinor_generators` and `spinor_genera_in_genus` (which
# enumerate_indefinite uses) gave 1 class instead of 2 in 13 of the 19162 indefinite
# genera of the 2026-10-22 sample; with the fix they agree with the patched Magma
# (magma_patches/diagonal-2.30-1.patch) on all of them.
#
# The fix is applied to the installed Hecke's own source text: the two functions are
# read from it, the PR's two changes are made, and the results are evaluated into
# Hecke.  If the source already has the fix, nothing is done; if it matches neither the
# buggy nor the fixed form, loading fails, since enumerating indefinite genera with the
# bug would silently store wrong class lists.

const HECKE_PR2611 = Ref(:not_checked)     # :applied or :already_in_hecke after __init__

function hecke_zgenus_source()
    for rel in (("QuadForm", "Quad", "ZGenus.jl"), ("QuadForm", "IntegerLattices", "ZGenus.jl"))
        path = joinpath(pkgdir(Hecke), "src", rel...)
        isfile(path) && return path
    end
    error("Hecke PR #2611: cannot find Hecke's ZGenus.jl under $(pkgdir(Hecke))")
end

# The top-level definition starting on the line that matches `start`, through the next
# line that is exactly "end".
function hecke_definition(lines::Vector{String}, start::Regex)
    i = findfirst(l -> occursin(start, l), lines)
    i === nothing && error("Hecke PR #2611: no definition matching $start in Hecke's ZGenus.jl")
    j = findnext(==("end"), lines, i)
    return join(lines[i:j], "\n")
end

function apply_hecke_pr2611!()
    lines = readlines(hecke_zgenus_source())
    a = hecke_definition(lines, r"^function automorphous_numbers\(g::ZZLocalGenus\)")
    b = hecke_definition(lines, r"^@attr Any function _automorphous_numbers\(G::ZZGenus\)")
    hunks = [(a, "  sort(I, by=x->valuation(x,2))", "  sort!(I, by=x->valuation(x,2))", "sort!(I"),
             (b, r"for p in P\n\s*push!\(gens_automorph, inj\[p\]\(QQ\(-1\)\)\)\n\s*end",
              "push!(gens_automorph, diagonal_map(QQ(-1)))", "diagonal_map(QQ(-1))")]
    applied = false
    for (code, buggy, fixed, marker) in hunks
        if occursin(marker, code)
            continue                                   # this part is already in Hecke
        end
        count(buggy, code) == 1 ||
            error("Hecke PR #2611: Hecke's source matches neither the buggy nor the fixed code; check the fix by hand")
        include_string(Hecke, replace(code, buggy => fixed), "Hecke PR #2611 (LMFDBLattices)")
        applied = true
    end
    HECKE_PR2611[] = applied ? :applied : :already_in_hecke
    return HECKE_PR2611[]
end
