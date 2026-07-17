Table name: `lat_genera`

This table stores lattices (free Z-modules with a nondegenerate symmetric inner product) up to local equivalence (also refered to as the genus of the lattice).

| Column | Type | Description |
| --- | --- | --- |
| [label](https://beta.lmfdb.org/Lattice/Labels) | text | The part of the label that is constant across a genus |
| [rank](https://beta.lmfdb.org/knowledge/show/lattice.dimension) | smallint | the rank of the lattice |
| nplus | smallint | the number of positive diagonal entries over R, so that a positive definite lattice has nplus equal to dimension |
| nminus | smallint | the number of negative diagonal entries over R, so that a positive definite lattice has nminus equal to 0 |
| [class_number](https://beta.lmfdb.org/knowledge/show/lattice.class_number) | smallint | size of the genus |
| [disc_abs](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | absolute value of determinant of Gram matrix |
| [disc_sign](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | smallint | sign of determinant of Gram matrix |
| [disc_radical](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | radical of determinant of Gram matrix |
| [disc_witt](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | bigint | Witt discriminant |
| [disc_geometric](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Geometric discriminant |
| [disc_quadratic](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Quadratic discriminant |
| [disc_half](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | float4 | Half discriminant |
| [disc_2adic_unit](https://beta.lmfdb.org/knowledge/show/lattice.discriminant) | smallint | 2-adic unit discriminant |
| bad_primes | integer[] | primes dividing the determinant |
| conway_symbol | text | the Conway symbol for the genus |
| dual_conway_symbol | text | the Conway symbol of the rescaled dual genus |
| [level](https://beta.lmfdb.org/knowledge/show/lattice.level) | bigint | level of lattice |
| is_even | boolean | whether the lattice is even |
| discriminant_group_invs | integer[] | Smith-style invariants for the discriminant group |
| discriminant_group_exponent | integer | Exponent of the discriminant group |
| discriminant_form | integer[] | Quadratic form on the discriminant group, as a symmetric matrix |
| adjacency_matrix | jsonb | A dictionary with primes as keys and flattened p-neighbor multi-edge adjaceny matrix as values |
| adjacency_polynomials | jsonb | A dictionary with primes as keys and factored characteristic polynomials of adjacency matrices as values (as a list of pairs (f,e) with f a list of integer coefficients and e an exponent) |
| mass | numeric[] | numerator and denominator of the mass: sum of 1/Aut(L) for L in the genus (only for definite genera) |
| theta_distinguishing_prec | smallint | The smallest nonnegative integer m with the property that if two lattices in the genus have theta series that are equal up to precision m then they are isometric |
| is_theta_distinguished | boolean | whether theta series distinguish lattices within this genus |
| hash_function | text | a string picking one of a set of possible hash functions.  The hash value stored for lattices within this genus is computed using this hash function |
| is_hash_distinguished | boolean | whether the chosen hash function distinguishes lattices within this genus |
| ambient_genus | text | if provided, the label of the genus from which this genus was enumerated: each lattice in this genus arises as the orthogonal complement of a vector in some lattice of the ambient genus (which lattice, and which vector, are recorded per-lattice in lat_lattices via ambient_lattice and orthogonal_complement).  The ambient lattice is not constant across the genus, so only the ambient genus is recorded here |
| scale | integer | the gcd of the entries of the Gram matrix for any lattice in this genus |
| rep | integer[] | a Gram matrix for a lattice within this genus.  Not stored if the lattices in this genus are all in lat_lattices |
| genus_hash | bigint | a hash value that can be easily computed from a lattice in this genus (should be collision-free) |