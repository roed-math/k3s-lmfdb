# Tests for oscar/fix_festi_veniani.jl on a tiny tree built in a temporary directory.
#   julia --project=oscar oscar/test_fix_festi_veniani.jl      (from lattices/)
# Exits 1 on a failure.

include(joinpath(@__DIR__, "fix_festi_veniani.jl"))

const HERE = dirname(@__DIR__)
failures = 0
check(cond, msg) = cond || (global failures += 1; println("FAIL: ", msg))

# The index straight from Oscar: [O(q_L) : image of O(L)].
function direct_index(g)
    L = integer_lattice(; gram = matrix(QQ, g))
    im, inj = image_in_Oq(L)
    return divexact(order(codomain(inj)), order(im)), order(im)
end

record(fmt, vals) = join([get(vals, k, "\\N") for k in fmt], "|")
flat(g) = "{" * join([g[i, j] for i in 1:size(g, 1) for j in 1:size(g, 2)], ", ") * "}"

function build_tree(root, genera)
    for f in ("genera_basic.format", "lat_basic.format")
        cp(joinpath(HERE, f), joinpath(root, f))
    end
    gfmt = read_format(joinpath(root, "genera_basic.format"))
    lfmt = read_format(joinpath(root, "lat_basic.format"))
    for (label, g, invs, fv, sidecar_im) in genera
        p = label_path(root, "genera_basic", label); mkpath(dirname(p))
        write(p, record(gfmt, Dict("label" => label, "rep" => flat(g), "discriminant_group_invs" => invs)) * "\n")
        p = label_path(root, "lattice_basic_data", "$label.1"); mkpath(dirname(p))
        write(p, record(lfmt, Dict("label" => "$label.1", "gram" => flat(g), "festi_veniani_index" => fv)) * "\n")
        if sidecar_im !== nothing
            p = label_path(root, "festi_veniani", label); mkpath(dirname(p))
            write(p, "$label.1|$sidecar_im\n")
        end
    end
    return lfmt
end

stored(root, lfmt, label) =
    split(first(eachline(label_path(root, "lattice_basic_data", "$label.1"))), "|")[findfirst(==("festi_veniani_index"), lfmt)]

# A1 + A1: Aut(L^#/L) = GL_2(F_2) has order 6, O(q_L) only swaps the two elements of
# norm 1/2, and O(L) induces that swap, so the index is 1 where the old code stored 3.
A1A1 = [2 0; 0 2]
# A2 + <6> and A1 + A3: small mixed discriminant forms.
A2_6 = [2 -1 0; -1 2 0; 0 0 6]
A1A3 = [2 0 0 0; 0 2 -1 0; 0 -1 2 -1; 0 0 -1 2]

cases = [("2.2.4.t1", A1A1, "{2, 2}"), ("3.3.18.t2", A2_6, "{3, 6}"), ("4.4.8.t3", A1A3, "{2, 4}")]
expected = Dict{String, ZZRingElem}()
images = Dict{String, ZZRingElem}()
for (label, g, invs) in cases
    expected[label], images[label] = direct_index(g)
end
check(expected["2.2.4.t1"] == 1, "A1 + A1 should have index 1, image_in_Oq gives $(expected["2.2.4.t1"])")
# |Aut| of Z/2+Z/2 = GL_2(F_2), Z/3+Z/6 = GL_2(F_3) x 1, Z/2+Z/4 (dihedral of order 8),
# Z/4+Z/12 = Aut(Z/4+Z/4) x Aut(Z/3), Z/30 (Euler phi)
for (invs, n) in [([2, 2], 6), ([3, 6], 48), ([2, 4], 8), ([4, 12], 96 * 2), ([30], 8)]
    check(abelian_aut_order(ZZ.(invs)) == n, "abelian_aut_order($invs) = $(abelian_aut_order(ZZ.(invs))), expected $n")
end

# 1. New fill output: sidecar with the image order, index \N.
mktempdir() do root
    lfmt = build_tree(root, [(l, g, i, "\\N", images[l]) for (l, g, i) in cases])
    check(main([root]) == 0, "main on new-style data")
    for (l, _, _) in cases
        check(stored(root, lfmt, l) == string(expected[l]), "new-style $l: $(stored(root, lfmt, l)) vs $(expected[l])")
    end
end

# 2. Legacy data: no sidecar, the column holds #Aut(L^#/L) / |image|; converted once.
mktempdir() do root
    legacy = Dict(l => string(divexact(abelian_aut_order(parse_list(i)), images[l])) for (l, _, i) in cases)
    lfmt = build_tree(root, [(l, g, i, legacy[l], nothing) for (l, g, i) in cases])
    check(stored(root, lfmt, "2.2.4.t1") == "3", "the old code's value for A1 + A1 is 3")
    check(main([root]) == 0 && stored(root, lfmt, "2.2.4.t1") == "3", "without --legacy, legacy data is left alone")
    for pass in 1:2                                  # the second pass must change nothing
        check(main(["--legacy", root]) == 0, "main --legacy, pass $pass")
        for (l, _, _) in cases
            check(stored(root, lfmt, l) == string(expected[l]), "legacy $l pass $pass: $(stored(root, lfmt, l)) vs $(expected[l])")
        end
    end
end

println(failures == 0 ? "test_fix_festi_veniani: all passed" : "test_fix_festi_veniani: $failures failures")
exit(failures == 0 ? 0 : 1)
