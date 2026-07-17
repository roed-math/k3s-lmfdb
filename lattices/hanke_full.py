from sage.all import *
import random
from collections import defaultdict
import itertools


def gcd_extended(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = gcd_extended(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return gcd, x, y

def find_min_x(num, rem):
    prod = 1
    for n in num:
        prod *= n

    result = 0
    for i in range(len(num)):
        prod_i = prod // num[i]
        _, inv_i, _ = gcd_extended(prod_i, num[i])
        result += rem[i] * prod_i * inv_i

    return result % prod

def algo3_8(L, a=2, at_prime=None):
    """a-saturate L (Hanke, arXiv:1208.2481, "Finding a Local a-saturated lattice").

    Globally by default; at_prime restricts it to a single prime.  L' is a-saturated when
    every Jordan scale is 0 or 1 -- i.e. L_p = L_0 + L_1 with L_i being a*p^i-modular.

    The algorithm is local and one-pass: Jordan-decompose L_p = (+)_i L_{i,p} into
    a*p^i-modular blocks and return (+)_i p^{-nu_i} L_{i,p} with nu_i = floor(i / 2).
    Scaling an a*p^i-modular block by p^{-floor(i/2)} makes it a*p^(i mod 2)-modular, which
    is why every scale lands in {0, 1} and why no iteration is needed.

    Restricted, the result is saturated at at_prime and its completions at every other prime
    q are unchanged (L'_q = L_q) -- what Sage's maximal_overlattice(p=...) promises, and
    what maximal_overlattice_2 relies on.
    """
    from sage.quadratic_forms.genera.normal_form import p_adic_normal_form

    IPM = L.inner_product_matrix()
    primes = ([Integer(at_prime)] if at_prime is not None
              else Integer(L.gram_matrix().determinant()).prime_factors())

    for p in primes:
        B = L.basis_matrix().change_ring(QQ)           # ambient coordinates
        G = L.gram_matrix()
        prec = Integer(G.determinant()).valuation(p) + 3
        # Jordan decomposition at p: p_adic_normal_form gives (D, T) with D = T*G*T^t block
        # diagonal over Z_p, the blocks a*p^i-modular -- exactly the L_i of the algorithm.
        D, T = p_adic_normal_form(G.change_ring(ZZ), p, precision=prec)
        # A block's scale i is the p-valuation of its rows.  Scaling an a*p^i-modular block
        # by p^{-floor(i/2)} leaves it a*p^(i mod 2)-modular, so "a-saturated" means every
        # scale is 0 or 1 -- and one pass suffices, no iteration.
        scales = []
        for i in range(D.nrows()):
            vals = [Integer(x).valuation(p) for x in D.row(i) if x != 0]
            scales.append(min(vals) if vals else 0)
        if all(s <= 1 for s in scales):
            continue                                   # already a-saturated at p
        S = diagonal_matrix(QQ, [QQ(1) / p**(s // 2) for s in scales])

        # Build L' directly as an OVERLATTICE rather than by transforming the basis.
        # Scaling a block by p^{-nu} enlarges the lattice inside the same Q_p-space, so
        # L'_p contains L_p, and no local-global gluing step is needed.
        #
        # Applying S*T as a global change of basis (the obvious move) is wrong: T is only
        # a p-adic transformation -- D = T*G*T^t holds mod p^prec -- so it corrupts every
        # other prime, and "works" only when det has a single prime factor.  Adjoining the
        # rows instead keeps the result honest at q != p: the denominators of S are powers
        # of p, hence q-units, so the new rows span the same Z_q-module as the rows of the
        # integral T, which already sits inside L_q.  Adding them back to L thus returns
        # L_q unchanged, while at p they span exactly the saturated L'_p.
        # T lives over Z_p.  Lift it through ZZ, never through QQ: change_ring(QQ) on a
        # p-adic matrix does *rational reconstruction*, so 1+2^2+2^4 (= 21) comes back as
        # -1/3 -- conjuring a denominator of 3 that then wrecks the lattice at q = 3.
        new = S * T.change_ring(ZZ).change_ring(QQ) * B
        L = IntegralLattice(IPM, Z_span_basis(list(B.rows()) + list(new.rows())))

    return L

def B(v, M, w=None):
    rng = M.base_ring()
    try:
        v = vector(rng, [rng(a) for a in v])
    except:
        pass
    if w != None:
        try:
            w = vector(rng, [rng(a) for a in w])
        except:
            pass
    if w is None:
        return (v * M * v.column())[0]
    else:
        return (v * M * w.column())[0]

def _to_field_vector(v, F):
    return vector(F, [F(x) for x in list(v)])

def B_field(v, M, w=None):
    F = M.base_ring()
    vF = _to_field_vector(v, F)
    if w is None:
        return (vF * M * vF.column())[0]
    wF = _to_field_vector(w, F)
    return (vF * M * wF.column())[0]


def radical_and_complement_rows_fp(M):
    """
    Returns (R_rows, Comp_rows) where rows are ambient row-vectors (as Matrix over M.base_ring()).
    """
    F = M.base_ring()
    n = M.nrows()

    # right kernel gives column vectors x with M*x = 0
    K = M.right_kernel()
    K_basis_cols = K.basis()   # list of column vectors (sage column vectors)
    # Radical rows: convert each column vector to a list of scalars
    R_rows = Matrix(F, [ list(col) for col in K_basis_cols ])  # may be 0 x n

    # Build a complement by taking standard basis columns and selecting those that extend span of radical
    base_cols = [ list(col) for col in K_basis_cols ]   # each entry is a length-n list
    comp_cols = []  # will store lists of scalars (each length n)

    for i in range(n):
        e = vector(F, [F(1) if j==i else F(0) for j in range(n)])
        e_list = list(e)
        # test matrix whose columns are base_cols + [e_list] + comp_cols
        test_cols_lists = base_cols + [e_list] + comp_cols
        base_cols_lists = base_cols + comp_cols
        if test_cols_lists:
            mat_test = Matrix(F, test_cols_lists).transpose()   # n x k
        else:
            mat_test = Matrix(F, 0, 0, [])
        if base_cols_lists:
            mat_base = Matrix(F, base_cols_lists).transpose()
        else:
            mat_base = Matrix(F, 0, 0, [])

        if mat_test.rank() > mat_base.rank():
            comp_cols.append(e_list)
        if len(comp_cols) + len(base_cols) >= n:
            break

    if comp_cols:
        Comp_rows = Matrix(F, comp_cols)   # rows = chosen complement vectors
    else:
        Comp_rows = Matrix(F, 0, n, [])

    return R_rows, Comp_rows

def restrict_gram(M, basis_rows):
    F = M.base_ring()
    if basis_rows.nrows() == 0:
        return Matrix(F, 0, 0, [])
    return basis_rows * M * basis_rows.transpose()

def _congruence_diagonalize(M):
    """Return (D, P) with D = P*M*P^t diagonal, over an odd-characteristic field.

    Symmetric Gaussian elimination: every row operation is mirrored by the same column
    operation (that is what keeps it a congruence, D = P M P^t, rather than a similarity),
    and accumulated in P.  Used by find_isotrop_fp to reduce isotropy to the diagonal case.
    """
    F = M.base_ring()
    n = M.nrows()
    A = copy(M)
    P = identity_matrix(F, n)
    for k in range(n):
        if A[k, k] == 0:
            piv = next((j for j in range(k + 1, n) if A[j, j] != 0), None)
            if piv is not None:
                A.swap_rows(k, piv); A.swap_columns(k, piv); P.swap_rows(k, piv)
            else:
                jj = next((j for j in range(k + 1, n) if A[k, j] != 0), None)
                if jj is None:
                    continue                      # row k is already zero from k on
                A.add_multiple_of_row(k, jj, 1)
                A.add_multiple_of_column(k, jj, 1)
                P.add_multiple_of_row(k, jj, 1)
        akk = A[k, k]
        if akk == 0:
            continue
        for j in range(k + 1, n):
            if A[j, k] != 0:
                f = -A[j, k] / akk
                A.add_multiple_of_row(j, k, f)
                A.add_multiple_of_column(j, k, f)
                P.add_multiple_of_row(j, k, f)
    return A, P


def find_isotrop_fp(M, max_trials=None):
    """A deterministic isotropic vector for q(v) = v*M*v^t over F_p, or None if the form is
    genuinely anisotropic.

    Complete by construction: it never reports anisotropic when an isotropic vector exists.
    The previous version sampled 200 random vectors and inferred anisotropy from failure,
    which over a large prime silently misses the isotropic vector that is there --
    maximal_overlattice_2(1009*A_2) returned det 3*1009^2 instead of 3, a non-maximal
    lattice with no exception raised.  max_trials is accepted and ignored (signature compat).
    """
    F = M.base_ring()
    n = M.nrows()
    V = VectorSpace(F, n)

    for i in range(n):
        if M[i, i] == 0:
            return V.gen(i)

    if F.characteristic() == 2:
        # Over F_2, q(v) = v*M*v^t = sum_i M_ii v_i is LINEAR (the cross terms 2*M_ij*v_i*v_j
        # vanish and v_i^2 = v_i).  Every diagonal entry is nonzero here (the zero case
        # returned above), so q(v) = sum_i v_i, and any weight-2 vector is isotropic.
        return V.gen(0) + V.gen(1) if n >= 2 else None

    # p odd: diagonalize, solve on the diagonal, map the solution back (x = y*P is isotropic
    # for M whenever y is isotropic for D = P*M*P^t, since x*M*x^t = y*D*y^t).
    D, P = _congruence_diagonalize(M)
    d = [D[i, i] for i in range(n)]
    for i in range(n):
        if d[i] == 0:
            return V.gen(i) * P
    # a coordinate 2-plane <d_i, d_j> is isotropic iff -d_j/d_i is a square (then
    # d_i*x^2 + d_j = 0 has a root); this settles every n = 2 case and many larger ones.
    for i in range(n):
        for j in range(i + 1, n):
            r = (-d[j]) / d[i]
            if r.is_square():
                y = [F(0)] * n
                y[i] = r.sqrt(); y[j] = F(1)
                return V(y) * P
    # n >= 3: the form is isotropic (any ternary form over F_p is), but possibly in no
    # coordinate 2-plane.  Solve the ternary d0*x^2 + d1*y^2 = -d2 by pigeonhole: {d0*x^2}
    # and {-d2 - d1*y^2} each hit (p+1)/2 values, so they meet.
    if n >= 3:
        a, b, c = d[0], d[1], d[2]
        squares = {}
        for x in F:
            squares.setdefault(a * x * x, x)
        for yv in F:
            t = -c - b * yv * yv
            if t in squares:
                y = [F(0)] * n
                y[0] = squares[t]; y[1] = yv; y[2] = F(1)
                return V(y) * P
    return None

def hyperbolic_fp(M):
    F = M.base_ring()
    n = M.nrows()
    V = VectorSpace(F, n)

    v = find_isotrop_fp(M)
    if v is None:
        return None

    for _ in range(3 * n + 10):
        w = V.random_element()
        if w == V.zero():
            continue
        val = B_field(v, M, w)
        if val != F(0):
            inv = val**(-1)
            w = inv * w
            if matrix(F, [list(v), list(w)]).rank() == 2:
                return [v, w]
    for e in V.basis():
        val = B_field(v, M, e)
        if val != F(0):
            w = (val**(-1)) * e
            if matrix(F, [list(v), list(w)]).rank() == 2:
                return [v, w]
    return None

def max_isotrop_fp(M, verbose=False):
    F = M.base_ring()
    n = M.nrows()

    R_rows, B_rows = radical_and_complement_rows_fp(M)
    radical_ambient = [ vector(F, list(r)) for r in R_rows.rows() ]

    if B_rows.nrows() == 0:
        return radical_ambient

    Gram_sub = restrict_gram(M, B_rows)
    planes = []
    while True:
        if B_rows.nrows() == 0 or Gram_sub.nrows() == 0:
            break

        out = hyperbolic_fp(Gram_sub)
        if out is None:
            break

        v_local, u_local = out   # local coordinates (length = Gram_sub.nrows())
        v_amb_row = (matrix(F, [list(v_local)]) * B_rows)[0]
        u_amb_row = (matrix(F, [list(u_local)]) * B_rows)[0]
        v_amb = vector(F, list(v_amb_row))
        u_amb = vector(F, list(u_amb_row))

        # use matrix * vector (no .column())
        w1 = M * v_amb
        w2 = M * u_amb
        W = Matrix(F, [ list(w1), list(w2) ])
        # The next working subspace is {v,u}^perp INTERSECTED with the current subspace
        # V_k = rowspan(B_rows), i.e. {v,u}^perp taken inside V_k.  W.right_kernel() alone is
        # {v,u}^perp in the FULL space, which wanders back into directions already split off:
        # on H (+) diag(1,1) over F_2 the kernel oscillates between the first two coordinates
        # and the last two and the loop never terminates.  Intersecting with V_k is what
        # makes the complement strictly shrink each pass.
        V_k = (F ** n).subspace(B_rows.rows())
        K = W.right_kernel().intersection(V_k)
        if K.dimension() == 0:
            B_rows = Matrix(F, 0, n, [])
            Gram_sub = Matrix(F, 0, 0, [])
            planes.append((v_amb, u_amb))
            break

        K_basis_cols = K.basis()
        B_rows = Matrix(F, [ list(col) for col in K_basis_cols ])
        Gram_sub = restrict_gram(M, B_rows)
        planes.append((v_amb, u_amb))
        if verbose:
            # print("Split off plane. Remaining dim (complement) =", B_rows.nrows())
            pass

    isotropic_list = radical_ambient + [ pair[0] for pair in planes ]
    return isotropic_list

def char2_max_isotrop(M):
    F = M.base_ring()
    n = M.nrows()
    V = VectorSpace(F, n)
    R_rows, B_rows = radical_and_complement_rows_fp(M)
    if B_rows.nrows() == 0:
        # only radical present
        return V.subspace([ vector(F, list(r)) for r in R_rows.rows() ])

    Gram_sub = restrict_gram(M, B_rows)
    rows = list(B_rows.rows())   # list of ambient row-vectors
    r = len(rows)

    # find a row whose self-pairing in Gram_sub is nonzero
    done = False
    for i, _row in enumerate(rows):
        if Gram_sub[i, i] != 0:
            done = True
            rows[i], rows[-1] = rows[-1], rows[i]
            break

    if done:
        # rebuild B_rows and Gram_sub after the swap
        B_rows = Matrix(F, [ list(rv) for rv in rows ])
        Gram_sub = restrict_gram(M, B_rows)
        rows = list(B_rows.rows())

        # try to make all but last have zero self-pairing by adding last row
        for i in range(len(rows) - 1):
            Gram_sub = restrict_gram(M, B_rows)
            if Gram_sub[i, i] != 0:
                rows[i] = rows[i] + rows[-1]
                B_rows = Matrix(F, [ list(rv) for rv in rows ])

    # collect those rows whose self-pairing (in Gram_sub) is zero
    Gram_sub = restrict_gram(M, B_rows)
    works = []
    for i, row in enumerate(B_rows.rows()):
        if Gram_sub[i, i] == 0:
            works.append(vector(F, list(row)))

    # build subspace: radical + the "works" rows (ambient vectors)
    ambient_rows = [ vector(F, list(r)) for r in R_rows.rows() ] + works
    return ambient_rows

def L_perp_mod2_basis(G, w):
    n = G.nrows()
    w = vector(ZZ, w)
    a = G * w  
    a = [x%2 for x in a]
    vecs = []
    bad = []
    for i in range(n):
        if a[i] == 0:
            zero = [0] * n
            zero[i] = 1
            vecs.append(zero)
        else:
            bad.append(i)
    curr = [0] * n
    for i in range(len(bad)):
        curr[bad[i]] = 1
        if i%2 == 0:
            curr[bad[0]] = 2
        vecs.append(curr.copy())
        curr[bad[0]] = 1

    return vecs

def Z_span_basis(gens):
    """Basis of the Z-module spanned by gens, via Hermite form.

    The generators are the ROWS: the previous version built matrix(QQ, n, len(gens), gens),
    which is transposed and only happens to have the right entry count -- p_neighbor_lattice
    already worked around it with an inline HNF rather than calling this.

    This is also the fast path.  Sage's overlattice()/span() compute the same span through
    the generic free_module machinery, which profiling showed to be ~95% of the cost of
    maximal_overlattice_2; doing it directly on integer matrices measured 18-23x faster
    with identical determinants.
    """
    Mg = matrix(QQ, [list(g) for g in gens])          # len(gens) x n, generators as rows
    den = lcm([c.denominator() for c in Mg.list()])
    H = (den * Mg).change_ring(ZZ).hermite_form(include_zero_rows=False)
    return H / den


def p_neighbor_lattice(L_in, w, p=2):
    G = L_in.gram_matrix()

    perp_basis = L_perp_mod2_basis(G, w)    

    gens = [vector(QQ, w) / p] + [vector(QQ, v) for v in perp_basis]
    # B = Z_span_basis(gens)
    mat = Matrix(gens)
    int_mat = mat.denominator() * mat
    int_mat = int_mat.change_ring(ZZ)
    hnf_int = int_mat.hermite_form(include_zero_rows=False)
    B = hnf_int / mat.denominator()
    Gprime = B * G * B.transpose()

    return IntegralLattice(Gprime)

def even_sublattice(L):
    G = L.gram_matrix()
    n = G.nrows()
    d = vector(ZZ, G.diagonal()) % 2
    if d.is_zero():
        return L
    pivot = []
    basis = []
    for i in range(n):
        if d[i] == 1:
            pivot.append(i)
        else:
            v = [0] * n
            v[i] = 1
            basis.append(v)
    v = [0]*n
    v[pivot[0]] = 2
    basis.append(v.copy())
    for i in range(1,len(pivot)):
        v[pivot[i]] = 1
        v[pivot[0]] = 3 - v[pivot[0]]
        basis.append(v.copy())
    B = Matrix(ZZ, basis)
    B = B * L.basis_matrix()
    # print(B)
    #return [IntegralLattice(B * L.gram_matrix() * B.transpose()), B]
    return L.sublattice(B)

# ------------------------
# Integer kernel basis for B v == 0 (mod 2), and finishing logic
# ------------------------

def integer_basis(B):
    B = Matrix(ZZ, B)   
    n = B.ncols()
    B2 = B.change_ring(GF(2))
    K = B2.right_kernel()
    basis_gf2 = K.basis()   
    lifted = []
    for vec in basis_gf2:
        lifted.append(vector(ZZ, [int(vec[j]) for j in range(n)]))

    two_std = [2 * vector(ZZ, [1 if i == j else 0 for j in range(n)]) for i in range(n)]

    Free = ZZ**n       
    Ssub = Free.submodule(lifted + two_std)
    zbasis = Ssub.basis()  

    return zbasis

def characteristic_vector(G):
    """A characteristic vector w: q(x) = b(x,w) mod 2 for every x, i.e. G*w = diag(G) mod 2.

    In characteristic 2 the cross terms 2*G_ij*x_i*x_j vanish and x_i^2 = x_i, so
    q(x) = sum_i G_ii*x_i mod 2 -- a LINEAR functional.  It is therefore b(x,w) for a
    single w, obtained by solving G*w = diag(G) over F_2.  Returns None when the system
    has no solution (b need not be onto mod 2 once det is even).
    """
    n = G.nrows()
    G2 = G.change_ring(GF(2))
    d = vector(GF(2), [G[i, i] for i in range(n)])
    try:
        w2 = G2.solve_right(d)
    except ValueError:
        return None
    return vector(ZZ, [Integer(int(x)) for x in w2])


def even_sublattice_at(L, p=2):
    """L_{a-even at p} = {x in L : H(x,x) in p^ord_p(2a)}, for a = Z.

    Hanke, Sec. 4 ("Even bilinear lattices") distinguishes the FULL even sublattice
    L_{a-even} = {x : H(x,x) in 2a} from the PARTIAL one at p, and notes
    L_{a-even} = intersection over p of L_{a-even at p}.  Over Z the only prime that can
    fail evenness is 2, so the two coincide and even_sublattice(L) is exactly this; the
    name is kept to match the algorithm's step 2.
    """
    return even_sublattice(L)


def residual_quadric_points(L, p=2):
    """Non-singular points of the residual a-quadric C_{p;a} at p, as lifts v_P in L.

    Hanke, Sec. 5:
      C_{p;a} := P({x in L/pL : Q(x) in a*p, x != 0})            (residual quadric)
      P is NONSINGULAR  <=>  grad Q(P) != 0 in (a/pa)^n
                        <=>  the linear form H(v_P, .) on L/pL is not identically zero
                             (Lemma "Hessian singularity criterion", since H(x,y) = grad Q(x).y)

    Yields v_P sheared so that Q(v_P) is in a*p^2, which the theorem requires ("v_P is any
    lift of a non-zero isotropic vector in L/p^2 L, i.e. Q(v_P) in p^2 a").  The shear is
    the proof's step 1: non-singularity gives a lift w of grad Q(v_P) with H(v_P, w) not in
    pa, so some lambda in p makes Q(v_P + lambda*w) land in a*p^2.
    """
    G = L.gram_matrix()
    n = G.nrows()
    Fp = GF(p)

    # Only CHARACTERISTIC vectors can give an a-even p-neighbour, so search those rather
    # than all p^n points of P(L/pL).  N(v) is a-even iff q vanishes on v^perp mod 2;
    # q(.) mod 2 is the linear form b(., w) (in characteristic 2 the cross terms drop and
    # x_i^2 = x_i), so it vanishes on v^perp exactly when b(., v) and b(., w) cut out the
    # same hyperplane, i.e. when G*v = diag(G) mod 2.  Those v form a coset of
    # ker(G mod 2): 2^corank candidates instead of 2^n.  For an odd L each is automatically
    # non-singular, since G*v = diag(G) != 0 mod p is exactly the theorem's condition.
    G2 = G.change_ring(Fp)
    w = characteristic_vector(G)
    if w is None:
        return                             # no characteristic vector => no even neighbour
    v0 = vector(Fp, [Fp(int(x)) for x in w])
    K = G2.right_kernel()                  # the solutions form the coset v0 + K
    for k in K:
        v = vector(ZZ, [Integer(int(x)) for x in (v0 + k)])
        if not any(x % p for x in v):
            continue                       # must be non-zero in L/pL
        grad = G * v
        if all(Integer(x) % p == 0 for x in grad):
            continue                       # singular: the theorem gives no neighbour here
        # shear (proof step 1): non-singularity gives a lift w of grad Q(v) with
        # H(v, w) not in pa, so some lambda in p lands Q(v + lambda*w) in a*p^2 --
        # i.e. H = 0 mod 2p^2, which for p=2 is the "isotropic mod 8" condition.
        w = vector(ZZ, [Integer(x) for x in grad])
        for lam in range(0, p * p, p):
            u = v + lam * w
            if (u * G * u.column())[0] % (2 * p * p) == 0:
                yield u
                break


def finish(L):
    """Step 2-3 of Hanke's "Construct a maximal a-valued quadratic lattice" (Sec. 4).

    Given L maximal in the Hessian BILINEAR space, the algorithm reads:

        1. find an a-maximal bilinear lattice L; let S be the primes where L is not
           a-even at p  (necessarily p | 2)
        2. for each p in S, use the p-neighbour theorem to check if some p-neighbour L'
           of L is a-even at p.  If so replace L by L'; if not, replace L by L_{a-even at p}
        3. the remaining L is maximal a-even, hence a maximal a-valued QUADRATIC lattice

    Over Z with a = Z, S is empty or {2}, so this is a single pass at p = 2.

    The p-neighbours are enumerated via the paper's Theorem ("p-neighbors via residual
    geometry"): they are in bijection with the NON-SINGULAR points of the residual quadric
    at p, so the candidate vectors are exactly those with Q(v) in a*p^2 -- H(v,v) = 0 mod 8
    for p = 2 -- that are non-singular.  An earlier version replaced this with a
    characteristic-vector criterion (solve G*w = diag(G) mod 2, then test q(w) = 0 mod 8).
    That is only valid when G is invertible mod 2, i.e. odd determinant: at det 50 the
    solve fails, no neighbour is found, and the even sublattice is returned (det 200) where
    an even neighbour exists at det 50.  Non-singularity, not characteristicity, is the
    condition the theorem actually requires.
    """
    if L == even_sublattice(L):
        return L                      # already a-even: nothing in S, done
    # p = 2 is the only prime of Z that can fail evenness
    for v in residual_quadric_points(L, p=2):
        Lp = p_neighbor_lattice(L, v)
        if Lp == even_sublattice(Lp):
            return Lp                 # an a-even 2-neighbour: keeps the determinant
    return even_sublattice_at(L, p=2)  # none exists: fall back to L_{a-even at 2}

# ------------------------
# Helpers for reducing rationals mod p and building p*M^{-1} over F_p
# ------------------------

def _Q_to_Fp_entry(r, p):
    r = QQ(r)
    if r == 0:
        return GF(p)(0)
    a = ZZ(r.numerator())
    b = ZZ(r.denominator())
    va = a.valuation(p) if a != 0 else +Infinity
    vb = b.valuation(p) if b != 0 else 0
    # cancel common p-powers
    t = min(va if va != +Infinity else va, vb)
    if t != +Infinity and t > 0:
        a //= (p**t)
        b //= (p**t)
        va -= t
        vb -= t
    # if denominator still has p, valuation(r) > 0, so 0 mod p
    if b.valuation(p) > 0:
        return GF(p)(0)
    return GF(p)(a % p) * GF(p)(Integer(b % p)).inverse()

def _matrix_Q_to_Fp(MQ, p):
    return matrix(GF(p), [[_Q_to_Fp_entry(MQ[i,j], p) for j in range(MQ.ncols())] for i in range(MQ.nrows())])


def is_overlattice(L1, L2):
    B1 = matrix(QQ, L1.basis())     
    B2 = matrix(QQ, L2.basis())

    coords = B2 * B1.inverse() 

    return all(c in ZZ for c in coords.list())

def install_fast_maximal_overlattice():
    """Substitute maximal_overlattice_2 for Sage's IntegralLattice.maximal_overlattice.

    Sage's representative() -- the census bottleneck (genus.py:623) -- spends most of its
    time in maximal_overlattice, calling it directly and again once per prime inside
    local_modification, so substituting the one primitive reaches both hot spots and speeds
    up representative() without reimplementing it.

    Correctness: a 184-case sweep of the patched representative() (rank 3-16, signatures
    definite and indefinite, determinants 1-196, checked against Sage's own criterion
    Genus(rep) == genus) passes 184/184 -- no wrong answers, no hangs.  That sweep covers
    only SMALL primes, though; three separate bugs were fixed to get the substitution
    trustworthy, and the last two were found only after the sweep:
      * algo3_8 lifted a p-adic matrix through QQ (rational reconstruction turned 21 into
        -1/3, corrupting other primes);
      * max_isotrop_fp looped forever on H (+) diag(1,1) over F_2 (it took a plane's
        complement in the full space instead of within the current subspace);
      * find_isotrop_fp searched for an isotropic vector by 200 RANDOM trials and inferred
        anisotropy from failure -- complete for small primes, but over a large prime it
        misses the vector that is there, so maximal_overlattice_2(1009*A_2) returned a
        non-maximal det 3*1009^2 (Sage: 3) with no exception.  It is now deterministic
        (diagonalize, solve the diagonal form exactly), so completeness no longer depends
        on the prime.
    The lesson the sweep taught: a small-determinant pass is not evidence of correctness at
    large primes.  If this is ever wired into the live pipeline, extend the regression set
    to large-prime and high-rank determinants first.

    Speed: a 40-genus benchmark (rank 8-24, det up to 32768) measured 1.64x overall, 1.30x
    to 2.52x per case, never slower -- biggest at low rank / high det, softening as rank
    climbs and more time moves into Sage internals the patch does not touch.

    It MONKEY-PATCHES a Sage class, changing behaviour for every caller in the process, so
    it is opt-in rather than done at import.  Returns the original method so a caller can
    restore it.
    """
    import sage.modules.free_quadratic_module_integer_symmetric as fqm
    cls = fqm.FreeQuadraticModule_integer_symmetric
    original = cls.maximal_overlattice

    def _fast(self, p=None):
        return maximal_overlattice_2(self, p=p, do_asserts=False)

    cls.maximal_overlattice = _fast
    return original


def maximal_overlattice_2(L_in, p=None, do_asserts=True):
    """Maximal INTEGRAL overlattice of L_in, with odd results allowed.

    p restricts the construction to a single prime: the result is maximal at p and
    unchanged at every other prime.  This mirrors Sage's maximal_overlattice(p=...) --
    whose only difference from the p=None case is likewise which primes it loops over --
    so that this can stand in for it at both call sites, including the
    maximal_overlattice(p=p) inside local_modification.

    Which notion this is matters, because there are two.  The isotropy test below is the
    one for the discriminant BILINEAR form -- adjoin a class when b(v,w) stays integral --
    so the result is maximal among *integral* overlattices and may be odd.  Sage's
    IntegralLattice.maximal_overlattice instead returns the maximal *EVEN* overlattice,
    which is a strictly stronger condition (q(v) in 2Z, not just b integral): at rank 12
    det 100 Sage stops at det 4 while this reaches det 1.  Neither is wrong; they answer
    different questions, and the determinants are not comparable.

    Note the pipeline needs the EVEN notion for now: Sage's representative() calls
    maximal_overlattice internally and assumes evenness, so this is not a drop-in for it.
    Getting the even notion out of this construction means changing the p = 2 test from
    "b(v,w) integral" (an F_2 condition) to "q(v) in 2Z" (a Z/4 condition, since for
    2-elementary D the values q(v) live in (1/2)Z and 2*q(v) mod 4 is what decides it).
    Odd p is unaffected: 2 is invertible there, so parity cannot change.
    """
    ogL = L_in
    L = L_in

    # Step 1: Z-saturate -- restricted to p when we were asked for a single prime, so the
    # completions at every other prime stay put (local_modification requires that).
    L_sat = algo3_8(L, at_prime=p)
    # Step 2: Work prime-by-prime on the dual to adjoin isotropic classes from D(L)
    M = L_sat.gram_matrix()
    Minv = M.inverse()   # over QQ
    detM = Integer(M.determinant())
    if p is None:
        ps = detM.prime_factors()
        if 2 not in ps:
            ps.append(2)          # 2 is always visited: finish() settles parity there
    else:
        ps = [Integer(p)]         # maximal at p only, other primes untouched

    # Adjoin one prime at a time, recomputing the Gram in between, rather than pooling
    # every prime's vectors into a single overlattice() call.  The F_p model below only
    # constrains pairs coming from the SAME p: it certifies b(v,w) integral for v,w in
    # one prime's isotropic subspace, but says nothing about b(v_2, v_5).  Those lifts
    # v*M^{-1} are arbitrary elements of L*, not of the p-primary component D_p, so the
    # usual "D_2 and D_5 are orthogonal" argument does not apply to them and mixing the
    # primes produced a non-integral Gram (every multi-prime genus tested failed:
    # D = (10,10) always, D = (2,2) never).  Handling each prime against the current
    # lattice is also correct because adjoining p-power-index classes leaves the other
    # primes' discriminant parts untouched.
    # Stay in matrix-land across the loop: carry the ambient basis B and Gram M, and build
    # an IntegralLattice once at the end.  Profiling showed ~95% of this function's time
    # was Sage's generic module layer (free_module.__init__, span, overlattice,
    # _element_constructor_), not the mathematics -- max_isotrop_fp itself is ~0.06s/call.
    # Spanning via Hermite form on integer matrices instead measured 18-23x faster with
    # identical determinants, but only if we do not rebuild a lattice object every pass.
    IPM = L_sat.inner_product_matrix()
    B = matrix(QQ, [list(b) for b in L_sat.basis()])   # ambient coordinates
    M = B * IPM * B.transpose()

    for q in ps:
        Minv = M.inverse()
        # Exponent of the discriminant group D = Z^n/M Z^n, for the p-primary projection
        # below.  Its elementary divisors are those of M, so take them straight off the
        # matrix rather than paying for discriminant_group() (0.713s of a 1.536s profile).
        m_exp = Integer(1)
        for d in M.change_ring(ZZ).elementary_divisors():
            if d:
                m_exp = m_exp.lcm(Integer(d))
        to_adjoin = []
        # Model the p-primary discriminant form on F_p^n by the Gram p*M^{-1} mod p:
        # writing x = v*M^{-1} in L*, we have p*q(x) = v*(p*M^{-1})*v^t, so v is
        # isotropic for this form exactly when q(x) is integral, and likewise
        # b(v,w) = 0 exactly when the inner product of the lifts is integral.  The
        # radical corresponds to vectors already in L.
        Mp_dual = _matrix_Q_to_Fp(q*Minv, q)

        # max_isotrop_fp returns a *totally* isotropic subspace: its radical part is
        # isotropic (b(r,.) = 0 forces b(r,r) = 0), each plane vector is isotropic, and
        # each successive plane is split off inside the previous ones' orthogonal
        # complement -- so the returned vectors are pairwise orthogonal.  This is what
        # makes the adjoined overlattice integral, and it holds in every characteristic,
        # including 2 (where b(v,v) = v*M*v^t = q(v) since the cross terms 2*M_ij*v_i*v_j
        # vanish, so a b-isotropic subspace is automatically q-isotropic).
        #
        # It replaces char2_max_isotrop for p = 2, which collected every row of zero
        # SELF-pairing without ever checking b(v,w) for v != w.  Those vectors are
        # individually isotropic but do not span an isotropic subspace, so adjoining them
        # produced a non-integral Gram (a concrete rank-8 det-100 genus gave b(v,w) =
        # -3/2) and overlattice() rejected it.  Note also that no per-vector "q integral"
        # filter belongs here: integrality is a property of the whole subgroup, and
        # filtering a subspace vector-by-vector is exactly what destroys it.
        iso_list = max_isotrop_fp(Mp_dual, verbose=False)  # list of row vectors over F_p

        # Project each lift into the p-primary component D_p before adjoining it.
        # p*M^{-1} is p-integral, so reducing it mod p is p-ADIC reduction: the model
        # certifies only that b(v,w) is p-integral, not that it is integral.  A raw lift
        # v*M^{-1} is an arbitrary element of L*, so at det = 2^2*5^2 two vectors from the
        # p=2 model can pair to k/5 -- 2-integral, isotropic in the F_2 model, and not
        # integral (39 of 105 pairs failed this way on a rank-8 det-100 genus).  Scaling
        # by the prime-to-p part of the discriminant group's exponent lands the class in
        # D_p, where the only denominators are powers of p, so p-integral now means
        # integral and the model's guarantee is the one we need.  The cofactor is a unit
        # on D_p, so this reparametrises the isotropic subgroup rather than shrinking it.
        Fp = GF(q)
        cofactor = m_exp // (q**m_exp.valuation(q))
        for v in iso_list:
            vZ = vector(ZZ, [int(Fp(x)) for x in v])  # 0..p-1 reps
            v_dual = cofactor * (vector(QQ, vZ) * Minv)   # in L*, now p-primary
            to_adjoin.append(v_dual * B)                  # ambient coordinates

        if to_adjoin:
            # Span basis + adjoined classes by Hermite form, rather than overlattice().
            B = Z_span_basis(list(B.rows()) + to_adjoin)
            M = B * IPM * B.transpose()

    L_sat = IntegralLattice(IPM, B)     # build the lattice object once, at the end
    # finish() settles parity, which is a purely 2-adic matter, so it must not run when we
    # were asked to work at an odd prime: local_modification calls maximal_overlattice(p=p)
    # on lattices it deliberately leaves odd, and even-ifying there changes the lattice at
    # 2 when only p was meant to move (at p=5 it returned the even sublattice, det 4, where
    # Sage gives det 1).
    saturated = L_sat                 # the saturate-and-adjoin result, before finish
    did_finish = (p is None or p == 2)
    if did_finish:
        L_sat = finish(saturated)
    if do_asserts:
        # The saturation and class-adjoining phase only ever ENLARGES the lattice, so up to
        # here the result is a genuine overlattice of the input.  Do NOT compare it against
        # Sage's maximal_overlattice: they compute different things.  Sage returns the
        # maximal EVEN overlattice, while the isotropy condition used above is the one for
        # the discriminant BILINEAR form (b(v,w) integral), which allows odd results -- at
        # rank 12 det 100 Sage stops at det 4 (even) whereas this reaches det 1 (odd).
        assert is_overlattice(saturated, ogL), "saturation is not an overlattice of the input"
        if did_finish:
            # finish() (Hanke Sec. 4) restores EVENNESS: it returns the maximal a-even
            # lattice in L's rational quadratic space, reached via an even 2-neighbour or,
            # when none exists, the index-2 even sublattice.  A 2-neighbour is not an
            # overlattice of L, and the even sublattice of an odd input is a SUBlattice --
            # so evenness, not overlattice-ness, is the invariant here.  (This is the
            # classical Z^6 no-neighbour case, which must return the even sublattice D_6.)
            assert L_sat.is_even(), "finish() did not return an even lattice"

    return L_sat

#L = IntegralLattice(Matrix(QQ, [[4,0,0,0,0], [0,4,0,0,0], [0,0,8,0,0], [0,0,0,4,0], [0,0,0,0,8]]))
#L_max = maximal_overlattice_2(L)
#print(L_max)
#print(L_max.gram_matrix())

"""
[1/2   0   0 1/2   0]
[  0 1/2   0 1/2   0]
[  0   0 1/4 1/2 1/4]
[  0   0   0   1   0]
[  0   0   0   0 1/2]

this is the sage code basis, the gram is the smae
"""
