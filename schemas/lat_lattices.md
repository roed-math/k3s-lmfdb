Table name: `lat_lattices`

This table stores lattices (free Z-modules with a nondegenerate symmetric inner product) up to isomorphism.

Label: `dimension.signature.determinant.genus_spec.tiebreaker` where

- `genus_spec` is

  - ommitted if determinant is 1 and signature is not a multiple of 8
  - 0 for the even lattice and 1 for the odd in the case that determinant is 1 and signature is a multiple of 8
  - otherwise, for each p with p^2|determinant, we give the concatenated dimensions of the p,p^2,p^3... blocks in the Jordan decomposition (using lower case letters and then upper case letters if these dimensions are larger than 9; we separate different primes with periods.  Finally, we encode the rest of the genus information (sign, scale, oddity etc) into a single integer as described below, and append it to the Jordan information.  This integer will be even if the genus is even and odd if the genus is odd.
  - Label is in the format r.s.d.j_1.j_2....j_k.x, where
    - r is the rank of the lattice
    - s is n_plus, the number of 1s in the diagonalization over R.
    - d is the absolute value of the determinant
    - If p_1, ... , p_k are the primes whose squares divide `2*d` `(p_i^2 | 2*d)`, then
    j_1,...,j_k are corresponding rank decompositions of their Jordan forms, omitting the first, encoded in base 62 (digits 0-9, then lowercase a-z then uppercase A-Z)
    For example, if the pairs of (valuation, rank) appearing in the decomposition are (3, 1), (4,10), (6,37), it will be encoded as 01a0B (the 0s come from the fact that the rank at valuations 2 and 5 are 0).
    - The last component of the label, x, is a hexadecimal string whose bits represent the local data.
    - Let n_2 be the number of non-zero blocks in the Jordan decomposition at 2.
    - The least n_2 bits specify the types (I or II) of the non-zero blocks at 2.
    - From these, once can deduce the compartments and trains in the local symbol at 2, let c, t be their numbers.
    - The next 3*c bits represent the oddities of the compartments, with every 3 bits giving an oddity mod 8.
    - The next t bits represent the signs of the trains.
    - Finally, for every other prime p dividing d, in increasing order, if there are n_p non-zero blocks in the Jordan decomposition at p, we add n_p bits representing the signs of these blocks.
  - For the tiebreaker, we use lexicographic sorting by canonical Gram matrix for definite lattices.  For indefinite lattices, a tiebreaker is only needed in rank 2 or 3 (in rank 3 spinor genera provide a complete invariant).

| Column | Type | Description |
| --- | --- | --- |
| [label](https://beta.lmfdb.org/Lattice/Labels) | text | We're changing the label; see above |
| [aut_size](https://beta.lmfdb.org/knowledge/show/lattice.group_order) | numeric | size of automorphism group |
| aut_label | text | label for the automorphism group as an abstract group |
| aut_group | text | string storing the automorphism group as a matrix group |
| [rank](https://beta.lmfdb.org/knowledge/show/lattice.dimension) | smallint | the rank of the lattice |
| nplus | smallint | the number of positive diagonal entries over R, so that a positive definite lattice has nplus equal to dimension |
| nminus | smallint | the number of negative diagonal entries over R, so that a positive definite lattice has nminus equal to 0 |
| [disc_abs](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | absolute value of determinant of Gram matrix |
| [disc_sign](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | smallint | sign of determinant of Gram matrix |
| [disc_radical](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | radical of determinant of Gram matrix |
| [disc_witt](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | Witt discriminant |
| [disc_geometric](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Geometric discriminant |
| [disc_quadratic](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Quadratic discriminant |
| [disc_half](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Half discriminant |
| [disc_2adic_unit](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | smallint | 2-adic unit discriminant |
| bad_primes | integer[] | primes dividing the determinant |
| [class_number](https://beta.lmfdb.org/knowledge/show/lattice.class_number) | smallint | size of the genus |
| [density](https://beta.lmfdb.org/knowledge/show/lattice.density) | numeric | density of the lattice centered sphere packing (only for definite lattices) |
| [center_density](https://beta.lmfdb.org/knowledge/show/lattice.density) | numeric | density of the lattice centered sphere packing (only for definite lattices) |
| [covering_norm](https://beta.lmfdb.org/knowledge/show/lattice.covering_radius) | integer[] | the square of the minimum real number so that balls of this radius around lattice points cover space (stored as numerator and denominator) |
| deep_hole_count | 
| [hermite](https://beta.lmfdb.org/knowledge/show/lattice.hermite_number) | numeric | Hermite number (only for definite lattices) |
| [kissing](https://beta.lmfdb.org/knowledge/show/lattice.kissing) | bigint | kissing number, the number of minimal vectors (only for definite lattices) |
| shortest | integer[] | orbit representatives for the set of shortest vectors under the action of the automorphism group |
| [level](https://beta.lmfdb.org/knowledge/show/lattice.level) | bigint | level of lattice |
| [minimum](https://beta.lmfdb.org/knowledge/show/lattice.minimal_vector) | integer | length of shortest vector (only for definite lattices) |
| name | text | a string like "A2 2E8", often null |
| [theta_series](https://beta.lmfdb.org/knowledge/show/lattice.theta) | numeric[] | a vector, counting the number of representations of n (odd) or 2n (even) |
| theta_prec | integer | Absolute precision of the theta series |
| [gram](https://beta.lmfdb.org/knowledge/show/lattice.gram) | smallint[] | A gram matrix for this isometry class |
| gram_others | numeric[] | A list of additional gram matrices.  This may include human-preferred gram matrices other than the canonical Gram matrix (for E8 for example), or Gram matrices whose entries are too large to fit in an integer |
| gram_is_canonical | boolean | whether the Gram matrix is canonical (null if not definite) |
| counter | integer | the position of this lattice within its genus (1-based), matching the sort order used to assign the label tiebreaker; ambient_lattice references in other genera point at this value |
| ambient_lattice | integer | the counter of the lattice, within the genus recorded as ambient_genus in lat_genera, in which this lattice arises as the orthogonal complement of the vector(s) in orthogonal_complement (null when the lattice was not produced by such a descent) |
| orthogonal_complement | integer[] | a vector or list of vectors so that this lattice is isometric to the orthogonal complement of their span, inside the lattice referenced by ambient_lattice.  When the enumeration works up a chain of genera toward a unimodular lattice, the vector is given in the first lattice along the chain for which a Gram matrix is provided |
| [canonical_gram] | integer[] | Canonical form for the Gram matrix; currently only available for definite lattices |
| [genus_label](https://beta.lmfdb.org/Lattice/Labels) | text | The part of the label that is constant across a genus |
| conway_symbol | text | the Conway symbol for the genus |
| pneighbors | jsonb | a dictionary with primes as keys and a list of labels as values (the p-neighbors) |
| discriminant_group_invs | integer[] | Smith-style invariants for the discriminant group |
| discriminant_group_exponent | integer | Exponent of the discriminant group |
| festi_veniani_index | numeric | the index of the lattice automorphism group within the automorphism group of the discriminant group |
| is_even | boolean | whether the lattice is even |
| dual_label | text | the label for the rescaled dual lattice (may be null if the discriminant is too large) |
| dual_theta_series | numeric[] | the theta series of the rescaled dual lattice |
| theta_prec | integer | Absolute precision of the dual theta series |
| dual_hermite | numeric | the Hermite number of the rescaled dual lattice (only for definite lattices) |
| dual_kissing | bigint | the kissing number of the rescaled dual lattice (only for definite lattices) |
| dual_density | numeric | the center density of the rescaled dual lattice (only for definite lattices) |
| dual_det | numeric | the determinant of the rescaled dual lattice |
| dual_conway_symbol | text | the Conway symbol for the rescaled dual lattice |
| is_universal | boolean | whether the quadratic form represents all positive integers (only for definite lattices) |
| is_even_universal | boolean | whether the quadratic form represents all positive even integers (only for definite lattices) |
| is_regular | boolean | whether the quadratic form represents all positive integers represented by its genus (only for definite lattices) |
| universality | smallint | the largest positive integer n so that every positive definite lattice of rank n embeds in this lattice (only for definite lattices) |
| even_universality | smallint | the largest positive integer n so that every even positive definite lattice of rank n embeds in this lattice (only for definite lattices) |
| regularity | smallint | the largest positive integer n so that every positive definite quadratic form of rank n that embeds in this genus embeds in this lattice (only for definite lattices) |
| is_indecomposable | boolean | whether the lattice is (orthogonally) indecomposable |
| is_additively_indecomposable | boolean | whether the lattice is additively indecomposable |
| orthogonal_factors | text[] | the orthogonal decomposition of the lattice (given as a duplicate-free list of lattice labels, sorted in reverse by multiplicity, then by label) |
| orthogonal_multiplicities | smallint[] | multiplicies of the lattices in the orthogonal decomposition (a list of integers of the same length as orthogonal_factors, in a corresponding order) |
| tensor_decompositions | jsonb | A list of lists of pairs.  The overall list contains different decompositions as a tensor product; the first entry of each pair is a label and the second a multiplicity. |
| is_tensor_product | boolean | whether this lattice has a nontrivial decomposition as a tensor product |
| root_sublattice | text | the name of the root sublattice (e.g. "I4 A2 2E8") |
| root_complement | text | the label for the orthogonal complement of the root sublattice |
| even_sublattice | text | the label for the sublattice generated by vectors of even norm |
| even_complement | text | the label for the orthogonal complement of the even sublattice |
| norm1_rank | smallint | the rank of the sublattice spanned by vectors of norm 1 |
| norm1_complement | text | the label for the complement of the norm 1 sublattice |
| successive_minima | integer[] | the sequence of successive minima, of length equal to the rank |
| scale | integer | the gcd of the entries of the Gram matrix |
| is_well_rounded | boolean | whether the minimal vectors span a finite index sublattice |
| is_minimal_vector_generated | boolean | whether the minimal vectors span the lattice |
| is_strongly_well_rounded | boolean | whether the minimal vectors contain a basis for the lattice |
| is_eutactic | boolean | whether the norm can be written as a positive combination of squared pairings with all of the minimal vectors |
| is_perfect | boolean | whether the space of real symmetric matrices is spanned by rank one matrices derived from the minimal vectors of the lattice |
| perfection_defect | smallint | the difference between n(n+1)/2 and the dimension of the space spanned by rank one matrices derived from the minimal vectors of the lattice |
| is_algebraic | boolean | Whether this lattice embeds in a number field, with pairing induced from the trace pairing |
| primitive_scaling | text | the label of the primitive lattice that is similar to this one (null if this lattice is already primitive) |
| t_design | smallint | the largest even integer t such that S is a spherical t-design (sum_{s in S} (x.s)^t = C * x.x^(t/2) for some C, which must be (min(L) #S)/n).  See Venkov's paper "reseaux et designs spheriques" |
| is_chiral | boolean | whether the automorphism group is contained in SO(n) |
| hash | bigint | the value of the genus-specific hash function on this lattice |