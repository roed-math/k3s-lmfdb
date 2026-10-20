import os
from functools import reduce

from sage.arith.misc import kronecker, prime_divisors, gcd # type: ignore
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors # type: ignore
from sage.interfaces.magma import magma # type: ignore
from sage.matrix.constructor import matrix # type: ignore
from sage.misc.functional import is_even, is_odd # type: ignore
from sage.misc.misc_c import prod # type: ignore
from sage.parallel.decorate import parallel # type: ignore
from sage.quadratic_forms.genera.genus import (
    GenusSymbol_global_ring,
    Genus_Symbol_p_adic_ring,
    LocalGenusSymbol,
) # type: ignore
from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring # type: ignore
from sage.rings.integer_ring import ZZ # type: ignore
from sage.rings.real_mpfr import RR # type: ignore
from hanke_full import install_fast_maximal_overlattice # type: ignore

def get_product(set_list):
    '''
    Returns all elements in the cartesian product of the sets in set_list

    :param set_list: list
    :return: list

    >>> get_product([[[1,2],[3,4,5]], [[6,7,8,9],[10]], [[11,12],[13,14,15]]])
    [[[1, 2], [6, 7, 8, 9], [11, 12]],
    [[1, 2], [6, 7, 8, 9], [13, 14, 15]],
    [[1, 2], [10], [11, 12]],
    [[1, 2], [10], [13, 14, 15]],
    [[3, 4, 5], [6, 7, 8, 9], [11, 12]],
    [[3, 4, 5], [6, 7, 8, 9], [13, 14, 15]],
    [[3, 4, 5], [10], [11, 12]],
    [[3, 4, 5], [10], [13, 14, 15]]]

    >>> get_product([])
    [[]]

    >>> get_product([[1]])
    [[1]]
    '''
    if len(set_list) == 0:
        return [[]]
    if len(set_list) == 1:
        return [[x] for x in set_list[0]]
    return [[x] + y for x in set_list[0] for y in get_product(set_list[1:])]

def get_dets(det_set, sign, num):
    '''
    Returns all possible sequences of length num of signs whose product is sign.
    These correspond to the quadratic character of the determinants of Jordan factors.

    :param det_set: list
    :param sign: int
    :param num: int

    :return: list

    >>> get_dets([1,-1],1,3)
    [[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]]
    '''
    dets = []
    b = len(det_set)
    num_opts = b**(num - 1)
    for n in range(num_opts):
        idxs = ZZ(n).digits(b, padto=num - 1)
        det = [det_set[idx] for idx in idxs]
        # sgn = sign*prod(det)**(-1)
        # Here we assume signs are 1 or -1
        sgn = sign*prod(det)
        det += [sgn]
        dets.append(det)
    return dets

def oddity(s):
    '''
    Compute the oddity for a list of quintuples that describe a 2-adic genus symbol
    '''
    odd = sum([t[-1] for t in s])
    k = len([t for t in s if t[0] % 2 == 1 and t[2] in [3,5]])
    odd += 4*k
    odd %= 8
    return odd

def expected_oddity(n_plus, n_minus, odd_symbols):
    '''
    Returns the oddity of the 2-adic genus symbol, given the signature and the p-adic genus symbols at the odd primes.
    This follows from the oddity condition.
    '''
    excess = sum([symb.excess() for symb in odd_symbols])
    return (excess + n_plus - n_minus) % 8

def all_local_genus_symbols(rank, det, p):
    '''
    Returns all possible p-adic genus symbols for a lattice of a given rank and determinant
    '''
    val_det, unit_det = ZZ(det).val_unit(p)
    symbs = []
    det_set = [1,-1]
    sign = unit_det.kronecker(p)
    val_det, unit_det = ZZ(det).val_unit(p)
    decs = WeightedIntegerVectors(val_det, range(1,val_det+1))
    ranks = [[rank-sum(dec)] + list(dec) for dec in decs if sum(dec) <= rank]
    for dec in ranks:
        scls = [i for i in range(len(dec)) if dec[i] > 0]
        rks = [dec[i] for i in scls]
        dets = get_dets(det_set, sign, len(scls))
        for det_list in dets:
            assert prod(det_list) == sign
            blocks = [[scls[j], rks[j], det_list[j]] for j in range(len(rks))]
            blocks.sort()
            symbs.append(blocks)
    return [Genus_Symbol_p_adic_ring(p,symb) for symb in symbs]

def listify_symbol(x):
    return reduce(lambda x,y : x+y, x.canonical_symbol())

def get_unique(symbs):
    '''
    Returns a list containing only one symbol from every equivalence class
    '''
    symb_dict = {}
    for i,s in enumerate(symbs):
        symb_dict[tuple(listify_symbol(s))] = i
    return [symbs[i] for i in symb_dict.values()]

def all_dyadic_genus_symbols(n_plus, n_minus, det, odd_symbols, only_even=True):
    '''
    Returns all the possible 2-adic genus symbols of a global genus symbol having a signature (n_plus, m_minus),
    determinant det and p-adic genus symbols odd_symbols at the odd primes p dividing the determinant
    '''
    rank = n_plus + n_minus
    p = 2
    val_det, unit_det = ZZ(det).val_unit(p)
    symbs = []
    det_set = [1,3]
    sign = unit_det.kronecker(p)
    decs = WeightedIntegerVectors(val_det, range(1,val_det+1))
    ranks = [[rank-sum(dec)] + list(dec) for dec in decs if sum(dec) <= rank]
    exp_oddity = expected_oddity(n_plus, n_minus, odd_symbols)
    for dec in ranks:
        scls = [i for i in range(len(dec)) if dec[i] > 0]
        rks = [dec[i] for i in scls]
        dets = get_dets([1,-1], sign, len(scls))
        for det_list in dets:
            assert prod(det_list) == sign
            blocks = [[scls[j], rks[j], det_set[ZZ(1 - det_list[j]) // 2]] for j in range(len(rks))]
            # We do not obtain all oddities.
            comps = []
            for block in blocks:
                oddities = []
                if (block[0] == 0) and only_even:
                    # if the scaling is trivial and we want an even lattice, we cannot have any oddity
                    oddities = [[0,0]]
                elif block[1] == 1:
                    if block[2] == 1: 
                        oddities = [[1,1],[1,7]]
                    else:
                        oddities = [[1,3],[1,5]]
                elif block[1] == 2:
                    if block[2] == 1:
                        oddities = [[0,0],[1,0],[1,2],[1,6]]
                    else:
                        oddities = [[0,0],[1,2],[1,4],[1,6]]
                elif block[1] % 2 == 0:
                    oddities = [[0,0],[1,0],[1,2],[1,4],[1,6]]
                else:
                    oddities = [[1,1],[1,3],[1,5],[1,7]]
                comps.append([block + odd for odd in oddities])
            # !! TODO - replace by something that generated only those with the correct oddity
            for cmps in get_product(comps):
                if (oddity(cmps) - exp_oddity) % 8 == 0:
                    symbs.append(cmps)
    gen_symbs = [Genus_Symbol_p_adic_ring(p,symb) for symb in symbs]
    if only_even:
        assert all([s.is_even() for s in gen_symbs])
    return get_unique(gen_symbs)

def all_genus_symbols(n_plus, n_minus, det, only_even=True):
    '''
    Returns all the genus symbols of even lattices with signature (n_plus, n_minus) and determinant det

    :param n_plus: int
    :param n_minus: int
    :param det: int
    :param only_even: bool
    :return: list

    >>> all_genus_symbols(16,0,1)
    [Genus of
    None
    Signature:  (16, 0)
    Genus symbol at 2:    1^16]

    >>> all_genus_symbols(16,0,1,only_even=False)
    [Genus of
    None
    Signature:  (16, 0)
    Genus symbol at 2:    1^16,
    Genus of
    None
    Signature:  (16, 0)
    Genus symbol at 2:    [1^16]_0]
    
    >>> len(all_genus_symbols(5,0,2*512))
    63

    >>> [len(all_genus_symbols(17,0,2*D)) for D in range(1,10)]
    [1, 1, 1, 3, 1, 2, 1, 3, 3]
    '''
    if is_odd(n_minus) and (det > 0):
        return []
    if is_even(n_minus) and (det < 0):
        return []
    rank = n_plus + n_minus
    primes = ZZ(2*det).prime_divisors()
    odd_primes = primes[1:]
    global_symbs = []
    odd_symbs = get_product([all_local_genus_symbols(rank, det, p) for p in odd_primes])
    for odd_symbols in odd_symbs:
        dyadic = all_dyadic_genus_symbols(n_plus, n_minus, det, odd_symbols, only_even=only_even)
        for symb2 in dyadic:
            g = GenusSymbol_global_ring((n_plus, n_minus), [symb2] + odd_symbols)
            global_symbs.append(g)
    return global_symbs

def encode_rank(rk):
    assert rk >= 0
    if (rk < 10):
        return str(rk)
    if (rk < 10 + 26):
        return chr(ord('a') + (rk - 10))
    assert rk < 10 + 26*2
    return chr(ord('A') + (rk - 10 - 26))

def create_genus_label(genus_sym):
    '''
    Returns a string which is the LMFDB label of the genus symbol.
    Label is in the format r.s.d.j_1.j_2....j_k.x, where 
    r is the rank of the lattices
    s is n_plus
    d is the absolute value of the determinant
    If p_1, ... , p_k are the primes whose squares divide d (p_i^2 | d), then
    j_1,...,j_k are corresponding rank decompositions of their Jordan forms, omitting the first, encoded in base 62
    (digits 0-9, then lowercase a-z then uppercase A-Z)
    For example, if the pairs of (valuation, rank) appearing in the decomposition are (3, 1), (4,10), (6,37), 
    it will be encoded as 01a0B.
    The last component of the label, x, is a hexadecimal string whose bits represent the local data.
    Let n_2 be the number of non-zero blocks in the Jordan decomposition at 2.
    The least n_2 bits specify the types (I or II) of the non-zero blocks at 2.
    From these, once can deduce the compartments and trains in the local symbol at 2, let c, t be their numbers.
    The next 3*c bits represent the oddities of the compartments, with every 3 bits giving an oddity mod 8.
    The next t bits represent the signs of the trains.
    Finally, for every other prime p dividing d, in increasing order, if there are n_p non-zero blocks in the 
    Jordan decomposition at p, we add n_p bits representing the signs of these blocks.

    :param genus_sym: GenusSymbol_global_ring
    :return: str 

    >>> q = QuadraticForm(ZZ,3, [1,0,1,1,0,3]).Hessian_matrix()
    >>> create_genus_label(Genus(q))
    '3.3.22.a6'

    >>> s = all_genus_symbols(5,0,2*61)[0]
    >>> create_genus_label(s)
    '5.5.122.66'

    >>> s = all_genus_symbols(8,0,2**5*3**4*5**3*7**2,only_even=False)[0]
    >>> create_genus_label(s)
    '8.8.15876000.00001.0001.001.01.a267'
    '''
    rk = genus_sym.rank()
    sig = genus_sym.signature() # the difference n_+ - n_-
    nplus = (rk + sig) // 2
    det = abs(genus_sym.determinant())
    primes = (2*det).prime_divisors()
    local_symbols = [genus_sym.local_symbol(p).symbol_tuple_list() for p in primes]
    local_symbols_filtered = [genus_sym.local_symbol(p).symbol_tuple_list() for p in primes if det.valuation(p) > 1]
    max_vals = [s[-1][0] for s in local_symbols_filtered]
    rks = [[0 for i in range(max_val)] for max_val in max_vals]
    for i, s in enumerate(local_symbols_filtered):
        for t in s:
            if t[0] > 0:
                rks[i][t[0]-1] = t[1]
    jordan_ranks = ["".join([encode_rank(r) for r in rk]) for rk in rks]

    # we start by encoding the local symbol at 2
    
    s2 = genus_sym.local_symbol(2)
    block_n = s2.number_of_blocks()
    train_ends = [x[-1] for x in s2.trains()]
    assert train_ends[-1] == block_n-1

    # The following can be made more efficient, but I really want to keep track of the single bits in the encoding for now,
    # for debbugging purposes
    bits = []
    # compartments are in correspondece with subsets of {0..block_n-1}
    comparts = reduce(lambda x,y:x+y, s2.compartments(),[])
    compart_symbol = sum([1 << e for e in comparts])
    # note that the compart_symbol will be even if and only if the lattice is even
    compart_nbits = block_n
    bits += ZZ(compart_symbol).digits(2, padto=compart_nbits) 
    
    oddities = [sum([s2.canonical_symbol()[t][4] for t in cmpart]) % 8 for cmpart in s2.compartments()]
    oddities_symbol = sum([o*(1 << (3*i)) for i,o in enumerate(oddities)])
    oddities_nbits = 3*len(s2.compartments())
    bits += ZZ(oddities_symbol).digits(2, padto=oddities_nbits)  
    
    signs = [s2.canonical_symbol()[train[0]][2] for train in s2.trains()]
    bits += [ZZ(1-s) // 2 for s in signs]

    assert primes[0] == 2
    for p in primes[1:]:
        s = genus_sym.local_symbol(p)
        signs = [t[2] for t in s.symbol_tuple_list()]        
        bits += [ZZ(1-s) // 2 for s in signs]
        
    local_data = sum([bits[i] << i for i in range(len(bits))])
    label = ".".join([str(x) for x in [rk,nplus,det] + jordan_ranks + [hex(local_data)[2:]]])
    return label

def decode_rank(c):
    '''
    Returns a number between 0 and 61, representing a rank of a component in the Jordan decomposition.
    The number is represented by a character c using 0-9a-zA-Z.
    
    :param c: str
    :return: int

    >>> decode_rank('a')
    10
    
    >>> decode_rank('A')
    36

    >>> decode_rank('9')
    9
    '''
    assert ((ord('0') <= ord(c) <= ord('9')) or (ord('a') <= ord(c) <= ord('z')) or (ord('A') <= ord(c) <= ord('Z')))
    if (ord('0') <= ord(c) <= ord('9')):
        return ord(c) - ord('0')
    if (ord('a') <= ord(c) <= ord('z')):
        return ord(c) - ord('a') + 10
    assert (ord('A') <= ord(c) <= ord('Z'))
    return ord(c) - ord('A') + 10 + 26

def build_compartments_and_trains(symbols, num_blocks_2, compart_bits):
    '''
    Returns two lists of lists. 
    The first list consists of compartments of the dyadic symbol.
    Here, each compartment is a maximal interval where all factors are of scaled type I.
    The second list consists of a list of trains of the dyadic symbol.
    Here, each train is a maximal interval having the property that for each pair of adjacent forms
    at least one is of scaled type I.

    :param symbols: list
    :param num_blocks_2: int
    :param compart_bits: list
    :return: tuple

    '''
    compartments = []
    trains = []
    in_compartment = False
    train_start = 0
    for i in range(num_blocks_2):
        if i > 0:
            scale_diff = symbols[0][i][0] - symbols[0][i-1][0]
            if (scale_diff > 2) and (i > train_start):
                trains.append([j for j in range(train_start,i)])
                train_start = i
            
        if compart_bits[i] == 1:
            if in_compartment:
                if (i > 0) and (scale_diff > 1):
                    compartments.append([j for j in range(compart_start, i)])
                    compart_start = i
            if not in_compartment:
                compart_start = i
                in_compartment = True
                if (i > train_start) and (scale_diff > 1):
                    trains.append([j for j in range(train_start,i)])
                    train_start = i
        else: # compart_bits[i] == 0
            if in_compartment:
                in_compartment = False
                compartments.append([j for j in range(compart_start, i+1)])
                if (i > train_start) and (scale_diff == 2):
                    trains.append([j for j in range(train_start,i)])
                    train_start = i
            else:
                if (i > train_start):
                    trains.append([j for j in range(train_start,i)])
                    train_start = i
                
    if in_compartment:
        compartments.append([j for j in range(compart_start, num_blocks_2)])

    trains.append([j for j in range(train_start, num_blocks_2)])
    return compartments, trains


def genus_symbol_from_label(label):
    '''
    Returns a genus symbol corresponding to an LMFDB label.

    :param label: str
    :return: GenusSymbol_global_ring

    >>> all_syms = all_genus_symbols(8,0,2**5*3**4*5**3*7**2,only_even=False)
    >>> all([genus_symbol_from_label(create_genus_label(s)) == s for s in all_syms])
    True
    '''
    split_label = label.split(".")
    rk = ZZ(split_label[0])
    sig = ZZ(split_label[1])
    n_plus = (rk + sig) // 2
    n_minus = (rk - sig) // 2
    det = ZZ(split_label[2])
    primes = (2*det).prime_divisors()
    symbols = []
    for i,p in enumerate(primes):
        rank_dec = [decode_rank(c) for c in split_label[3+i]]
        rank_dec = [rk - sum(rank_dec)] + rank_dec
        symbols_p = [[i, rank_dec[i], 1] for i in range(len(rank_dec)) if rank_dec[i] > 0]
        symbols.append(symbols_p)

    # initializing the symbol at 2
    symbols[0] = [s + [0,0] for s in symbols[0]]
    num_blocks_2 = len(symbols[0])
    
    local_data = ZZ(split_label[3 + len(primes)], base=16)

    compart_bits = (local_data % (1 << num_blocks_2)).digits(2, padto=num_blocks_2)

    for i in range(num_blocks_2):
        symbols[0][i][3] = compart_bits[i]

    compartments, trains = build_compartments_and_trains(symbols, num_blocks_2, compart_bits)
            
    local_data >>= num_blocks_2
    
    oddity_bits = (local_data % (1 << (3*len(compartments)))).digits(2, padto=3*len(compartments))

    oddity_nums = [oddity_bits[3*i:3*(i+1)] for i in range(len(compartments))]
    oddities = [sum([(b[i] << i) for i in range(3)]) for b in oddity_nums]

    # updating oddities
    
    for j,c in enumerate(compartments):
        symbols[0][c[0]][4] = oddities[j]

    local_data >>= 3 * len(compartments)

    # train_signs = [i for i in range(len(trains)) if (local_data.bits()[i+offset] == 1)]

    train_sign_bits = (local_data % (1 << len(trains))).digits(2, padto=len(trains))

    # updating signs at 2
    
    for t,train in enumerate(trains):
        sign = (-1)**(train_sign_bits[t])
        symbols[0][train[0]][2] = 3 if sign == -1 else 1

    local_data >>= len(trains)

    for j in range(1, len(symbols)):
        n = len(symbols[j])
        sign_bits = (local_data % (1 << n)).digits(2, padto=n)
        for i in range(n):
            symbols[j][i][2] = (-1)**(sign_bits[i])
        local_data >>= n

    local_symbols = [Genus_Symbol_p_adic_ring(p, symbols[j]) for j,p in enumerate(primes)]
    return GenusSymbol_global_ring((n_plus, n_minus), local_symbols)

# This is based on the __repr__ method from local symbols
# !! TODO !! - should refactor the code there and then use it

def conway_symbol_dyadic(local_symbol):
    CS = local_symbol.canonical_symbol()
    trains = local_symbol.trains()
    comps = local_symbol.compartments()

    # removing the 0-index in case we omit it

    if CS[0][0] == 0:
        if trains[0] == [0]:
            trains = trains[1:]
        else:
            trains = [trains[0][1:]] + trains[1:]

        if len(comps) > 0:
            if comps[0] == [0]:
                comps = comps[1:]
            elif comps[0][0] == 0:
                comps = [comps[0][1:]] + comps[1:]

    CS_string = ""
    for train in trains:
        # mark the beginning of a train with a colon
        CS_string += " :"
        # collect the indices where compartments begin and end
        compartment_begins = []
        compartment_ends = []
        for comp in comps:
            compartment_begins.append(comp[0])
            compartment_ends.append(comp[-1])

        for block_index in train:
            if block_index in compartment_begins:
                # mark the beginning of this compartment with [
                CS_string += "["
            block = CS[block_index]
            block_string = "%s^{%s}" % (p**block[0], block[2] * block[1])
            CS_string += block_string
            if block_index in compartment_ends:
                # close this compartment with ]
                CS_string += "]"
                # the oddity belongs to the compartment
                # and is saved in its first block
                i = compartment_ends.index(block_index)
                compartment_start = compartment_begins[i]
                oddity = CS[compartment_start][4]
                CS_string += "_%s" % oddity
        # remove the first colon
        CS_string = CS_string[2:]
        # remove some unnecessary whitespace
        CS_string = CS_string.replace("  :", ":")
    return CS_string

# Here we do not want the trains and compartments, for the global symbol
def conway_symbol_local_part(local_symbol):
    '''
    Returns the part of the Conway symbol of the global genus symbol corresponding to the local symbol

    :param local_symbol: Genus_Symbol_p_adic_ring
    :return: str

    >>> s2A = Genus_Symbol_p_adic_ring(2, [[0,2,1,0,0],[1,2,1,1,2],[3,2,1,1,0]])
    >>> s2B = Genus_Symbol_p_adic_ring(2, [[0,2,1,0,0],[1,2,1,1,0],[3,2,1,1,2]])
    >>> conway_symbol_local_part(s2A)
    '2^{2}_{2}8^{2}_{I}'
    >>> conway_symbol_local_part(s2B)
    '2^{2}_{I}8^{2}_{2}'
    '''
    p = local_symbol.prime()
    CS_string = ""
    symbols = local_symbol.canonical_symbol()[:]
    if (p == 2) and (symbols[0][3] == 0):
        oddities = [i for i,s in enumerate(symbols) if s[4] != 0]
        zero_odd = len([s for s in symbols if (s[3] == 1) and (s[4] == 0)])
        if (len(oddities) > 0) and (zero_odd == 0):
            last_idx = oddities[-1]
            # we won't display the oddity when it can be inferred from the oddity formula
            symbols[last_idx][4] = 0

    symbols = [s for s in symbols if s[0] > 0]
    for s in symbols:
        CS_string += "%s^{%s}" % (p**s[0], s[2] * s[1])
        if (p == 2) and ((s[1] % 2 == 0) or (s[1] == 1)):
            if s[3] == 0:
                CS_string += "_{II}"
            else:
                if s[4] == 0:
                    CS_string += "_{I}"
                else:
                    CS_string += "_{" + str(s[4]) + "}"
    return CS_string

def conway_symbol(genus_symbol):
    '''
    Returns the conway symbol string of the genus symbol.
    '''
    roman_numeral = "II" if genus_symbol.is_even() else  "I"
    r,s = genus_symbol.signature_pair()
    prefix = roman_numeral + "_{" + str(r) + "," + str(s) + "}"
    local_part = "".join([conway_symbol_local_part(s) for s in genus_symbol.local_symbols()])
    if len(local_part) > 0:
        return prefix + "(" + local_part + ")"
    return prefix

def create_genus_entry(genus_symbol):
    '''
    Returns a dictionary with fields corresponding to the schema,
    and values populated by the genus symbol.
    '''
    table_row = {}
    table_row['label'] = create_genus_label(genus_symbol)
    table_row['rank'] = rank = genus_symbol.rank()
    table_row['nplus'], table_row['nminus'] = genus_symbol.signature_pair()
    disc = genus_symbol.determinant()
    table_row['disc_abs'] = abs(disc)
    table_row['disc_sign'] = disc.sign()
    table_row['bad_primes'] = disc.prime_divisors()
    table_row['disc_radical'] = prod(table_row['bad_primes'])
    table_row['disc_witt'] = (-1)**(rank*(rank-1)//2) * disc
    table_row['disc_geometric'] = RR(abs(disc)).sqrt()
    table_row['disc_quadratic'] = RR(disc) / 2**rank
    table_row['disc_half'] = table_row['disc_quadratic']
    if genus_symbol.is_even() and rank % 2 == 1:
        table_row['disc_half'] /= 2
    table_row['disc_2adic_unit'] = disc.val_unit(2)[1] % 8
    table_row['conway_symbol'] = conway_symbol(genus_symbol)
    table_row['level'] = level = genus_symbol.level()
    table_row['is_even'] = genus_symbol.is_even()
    disc_form = genus_symbol.discriminant_form() # Timing: this was slow
    invs = disc_form.invariants()
    table_row['discriminant_group_invs'] = invs
    # the trivial discriminant group (e.g. unimodular lattices) has no invariants;
    # its exponent is 1
    table_row['discriminant_group_exponent'] = invs[-1] if invs else 1
    disc_q = disc_form.gram_matrix_quadratic()
    den = disc_q.denominator()
    table_row['discriminant_form'] = (den*disc_q).list()
    if genus_symbol.signature() == genus_symbol.rank():
        mass = genus_symbol.mass() # Timing: this was slow
        table_row['mass'] = [mass.numerator(), mass.denominator()]
    else:
        table_row['mass'] = None
    # We store a representative for use by the magma code
    rep = genus_symbol.representative()
    table_row['rep'] = rep.list() # Timing: this was slow
    table_row['scale'] = gcd(table_row['rep'])
    dual = rep**(-1)
    dual = (dual.denominator() * dual).change_ring(ZZ)
    ddet = dual.det()
    dual_local_symbols = [LocalGenusSymbol(dual, p) for p in (2*ddet).prime_divisors()]
    dual_symbol = GenusSymbol_global_ring(genus_symbol.signature_pair(), dual_local_symbols, representative=dual)
    table_row['dual_conway_symbol'] = conway_symbol(dual_symbol)
    # Every trace_bound currently in the LMFDB of Gamma0(N) or Gamma(N,chi) for chi a quadratic character has trace bound less than 150.

    # sage: list(db.mf_newspaces.search({"trace_bound":{"$gte":100}}, ["label", "trace_bound"]))
    # [{'label': '3912.1.cp', 'trace_bound': 201}, {'label': '3600.2.x', 'trace_bound': 101}, {'label': '6336.2.m', 'trace_bound': 143}, {'label': '7056.2.k', 'trace_bound': 121}, {'label': '9248.2.a', 'trace_bound': 145}]
    # sage: list(db.mf_gamma1.search({"trace_bound":{"$gte":100}}, ["label", "trace_bound"]))
    # [{'label': '2997.1', 'trace_bound': 103}, {'label': '3423.1', 'trace_bound': 147}, {'label': '3549.1', 'trace_bound': 144}, {'label': '3584.1', 'trace_bound': 109}, {'label': '1792.2', 'trace_bound': 193}, {'label': '2560.2', 'trace_bound': 108}, {'label': '2570.2', 'trace_bound': 130}, {'label': '2816.2', 'trace_bound': 129}, {'label': '3072.2', 'trace_bound': 233}, {'label': '3584.2', 'trace_bound': 193}, {'label': '3840.2', 'trace_bound': 169}, {'label': '1792.3', 'trace_bound': 193}, {'label': '2560.3', 'trace_bound': 108}, {'label': '2816.3', 'trace_bound': 129}, {'label': '3072.3', 'trace_bound': 233}, {'label': '3584.3', 'trace_bound': 193}, {'label': '3840.3', 'trace_bound': 169}, {'label': '1792.4', 'trace_bound': 193}]
    return table_row

def create_lattice_entry(gram):
    # Basic invs (rank, nplus, det, disc, class_number, discriminant_group_invs, is_even)
    # Gram (gram, canonical_gram)
    # Labels (genus_label, label, name, conway_symbol)
    # Dual (dual_label, dual_theta_series, dual_hermite, dual_kissing, dual_density, dual_det, dual_conway)
    # Automorphism group (aut_size, aut_label, festi_veniani_index)
    # Packing info (density, hermite, kissing, minimum) -> Density, HermiteNumber, KissingNumber, Minimum
    # Representing integers (level, theta_series) ThetaSeries (TODO: figure out precision.  Max(trace_bound, 100, minimum+4) or something; going up to trace bound should mean we distinguish between lattices in the genus that have different theta series.  Compute trace bound in half integral weight using Shimura isomorphism)
    # Neighbors (pneighbors)
    pass

def fill_genus_entry(label, hecke_primes=[2]):
    '''
    The list of hecke_primes determines for which primes we compute
    the Hecke matrices and Hecke polynomials
    '''
    # These operations are time consuming
    genus_symbol = genus_symbol_from_label(label)
    # In high rank, the following call should be replaced by Noam's code
    reps = genus_symbol.representatives()
    

    genus_row['class_number'] = h = len(reps)

    lat = reps[0]
    hecke_mats = {}
    hecke_polys = {}

    if genus_symbol.signature() == genus_symbol.rank():
        for p in hecke_primes:
            hecke_mag = magma.AdjacencyMatrix(magma.Genus(magma.LatticeWithGram(lat)),p)
            hecke_mats[p] = [ZZ(x) for x in list(magma.Eltseq(hecke_mag))]
            hecke_mat = matrix(ZZ, h, h, hecke_mats[p])
            hecke_poly_fac = list(hecke_mat.charpoly().factor())
            hecke_polys[p] = [[fa[0].list(), fa[1]] for fa in hecke_poly_fac]
        genus_row['adjacency_matrix'] = hecke_mats
        genus_row['adjacency_polynomials'] = hecke_polys
    else:
        genus_row['adjacency_matrix'] = None
        genus_row['adjacency_polynomials'] = None

    return genus_row, lattice_rows

# The fields here are copied from the schema
COL_TYPE_LATTICE_GENUS = {'det': 'bigint',
 'disc': 'bigint',
 'id': 'bigint',
 'level': 'bigint',
 'rank': 'smallint',
 'nplus': 'smallint',
 'is_even': 'boolean',
 'discriminant_form': 'integer[]',
 'discriminant_group_invs': 'integer[]',
 'conway_symbol': 'text',
 'label': 'text',
 'rep': 'integer[]',
 'mass': 'numeric[]',
 'class_number': 'smallint',
 'adjacency_matrix': 'jsonb',
 'adjacency_polynomials': 'jsonb',
 'theta_prec' : 'smallint',
 'dual_conway_symbol' : 'text'}

COL_TYPE_LATTICE = {'det_abs': 'bigint',
 'det_sign': 'smallint',
 'det_radical': 'bigint',
 'disc': 'bigint',
 'dual_det': 'numeric',
 'dual_kissing': 'bigint',
 'id': 'bigint',
 'kissing': 'bigint',
 'level': 'bigint',
 'festi_veniani_index': 'numeric',
 'minimum': 'integer',
 'class_number': 'smallint',
 'rank': 'smallint',
 'nplus': 'smallint',
 'is_even': 'boolean',
 'discriminant_group_invs': 'integer[]',
 'gram': 'integer[]',
 'gram_is_canonical': 'boolean',
 'gram_others': 'integer[]',
 'pneighbors': 'jsonb',
 'aut_size': 'numeric',
 'density': 'numeric',
 'dual_density': 'numeric',
 'dual_hermite': 'numeric',
 'hermite': 'numeric',
 'dual_theta_series': 'numeric[]',
 'theta_series': 'numeric[]',
 'aut_label': 'text',
 'conway_symbol': 'text',
 'dual_conway': 'text',
 'dual_label': 'text',
 'genus_label': 'text',
 'label': 'text',
 'aut_group': 'text',
 'theta_prec': 'smallint',
 'dual_theta_prec': 'smallint',
 'name': 'text',
 'is_indecomposable': 'boolean',
 'is_additively_indecomposable': 'boolean',
 'orthogonal_factors': 'text[]',
 'orthogonal_multiplicities': 'smallint[]',
 'tensor_decompositions': 'jsonb',
 'is_tensor_product': 'boolean',
 'root_sublattice': 'text',
 'root_complement': 'text',
 'even_sublattice': 'text',
 'even_complement': 'text',
 'norm1_sublattice': 'text',
 'norm1_complement': 'text',
 'Zn_complement': 'text',
 'successive_minima': 'integer[]',
 'shortest': 'integer[]'}

FIELDS_LATTICE_GENUS = ["genera_basic.format", "genera_advanced.format"]
FIELDS_LATTICE = ["lat_basic.format", "lat_advanced.format"]

def write_header_to_file(fname, sep = "|", col_type=COL_TYPE_LATTICE_GENUS, fields=FIELDS_LATTICE_GENUS):
    # we want to have a well defined order, matching the entries
    #fields = open("genera_basic.format").read().split("|")
    #fields += open("genera_advanced.format").read().split("|")
    fields = reduce(lambda x,y: x+y, [open(f).read().split("|") for f in fields])

    header_lines = [sep.join(fields), sep.join([col_type[k] for k in fields]), "\n"]
    
    header = "\n".join(header_lines)
    
    f = open(fname, "w")
    f.write(header)
    f.close()

def value_to_postgres(val):
    if type(val) == list:
        return str(val).replace('[', '{').replace(']','}')
    if type(val) == tuple:
        return value_to_postgres(list(val))
    if type(val) == str:
        # return '"' + val + '"'
        return val
    if type(val) == dict:
        d = {}
        for k in val.keys():
            d[str(k)] = val[k]
        return str(d)
    if type(val) == bool:
        return "T" if val else "F"
    return str(val)

def write_entries_to_file(entries, fname, sep = "|", col_type=COL_TYPE_LATTICE_GENUS):
    # we want to have a well defined order, matching the entries
    # fields = sorted(list(col_type.keys()))
    with open("genera_basic.format") as f:
        fields = f.read().split("|")
    lines = [sep.join([value_to_postgres(entry[k]) for k in fields]) for entry in entries]
    output = "\n".join(lines)
    f = open(fname, "a")
    num = f.write(output)
    if (num > 0):
        f.write("\n")
    f.close()

def write_header_and_entries(entries, fname, sep = "|", col_type=COL_TYPE_LATTICE_GENUS):
    write_header_to_file(fname, sep = sep, col_type=col_type)
    write_entries_to_file(entries, fname, sep = sep, col_type=col_type)

def write_all_of_sig_up_to_det(n_plus, n_minus, det):
    '''
    Create data file with all genera of a certain signature with determinant up to det.
    '''
    if not os.path.exists("data"):
        os.makedirs("data")
    fname = "data/genera_signature_%s_%s_%s.tbl" % (n_plus, n_minus, det)
    write_header_to_file(fname)
    sgn = 1 if is_even(n_minus) else -1;
    for d in range(1, det+1):
        print("determinant = %s" % (sgn*d))
        syms = all_genus_symbols(n_plus, n_minus, sgn*d, only_even=False)
        entries = [create_genus_entry(s) for s in syms]
        write_entries_to_file(entries, fname)

def write_all_of_sig_between(n_plus, n_minus, lb_det, ub_det):
    '''
    Create data file with all genera of a certain signature with determinant between lb_det and ub_det
    '''
    if not os.path.exists("data"):
        os.makedirs("data")
    fname = "data/genera_signature_%s_%s_%s_%s.tbl" % (n_plus, n_minus, lb_det, ub_det)
    write_header_to_file(fname)
    sgn = 1 if is_even(n_minus) else -1;
    for d in range(lb_det, ub_det+1):
        # print("determinant = %s" % (sgn*d))
        syms = all_genus_symbols(n_plus, n_minus, sgn*d, only_even=False)
        entries = [create_genus_entry(s) for s in syms]
        write_entries_to_file(entries, fname)

@parallel(ncpus=128) # !!! change this according to platform
def write_all_of_sig_between_genera_basic(n_plus, n_minus, lb_det, ub_det):
    '''
    Create data files with all genera of a certain signature with determinant between lb_det and ub_det,
    one file for each genus
    '''
    # Speeds up representative() (genus.py:create_genus_entry, the listing stage's main
    # cost) ~1.64x by substituting maximal_overlattice_2 for Sage's own maximal_overlattice.
    # @parallel forks a fresh process per call, so this must be installed here rather than
    # at module import: install_fast_maximal_overlattice monkey-patches a Sage class for
    # every caller in the process, and genus.py is also imported by code that never wants
    # that (e.g. testGenus.py). See hanke_full.install_fast_maximal_overlattice's docstring
    # for the correctness sweep this relies on.
    install_fast_maximal_overlattice()
    folder_name = "genera_basic/%d/%d" % (n_plus + n_minus, n_plus)   # rank/nplus, matching LabelPath
    os.makedirs(folder_name, exist_ok=True)
    sgn = 1 if is_even(n_minus) else -1;
    for d in range(lb_det, ub_det+1):
        syms = all_genus_symbols(n_plus, n_minus, sgn*d, only_even=False)
        entries = [create_genus_entry(s) for s in syms]
        for genus in entries:
            fname = folder_name + "/%s" % genus['label']
            if os.path.exists(fname):
                os.remove(fname)
            write_entries_to_file([genus], fname)


def write_all_up_to_det(rank, det):
    '''
    Create data files with all genera of a certain rank with determinant up to det.
    '''
    for n_minus in range(rank // 2 + 1):
        n_plus = rank - n_minus
        print("signature = (%s,%s)" %(n_plus, n_minus))
        write_all_of_sig_up_to_det(n_plus, n_minus, det)

def determinant_condition(sp, det):
    p = sp.prime()
    sgn = prod([x[2] for x in sp.symbol_tuple_list()])
    _, a = det.val_unit(p)
    return kronecker(sgn,p) == a.kronecker(p)

def oddity_condition(genus_symbol):
    s = genus_symbol
    s2 = s.local_symbol(2)
    oddity_s = oddity(s2.symbol_tuple_list())
    excess = sum([sp.excess() for sp in s.local_symbols() if sp.prime() != 2])
    return (s.signature() + excess - oddity_s) % 8 == 0

def jordan_condition(genus_symbol):
    s2 = genus_symbol.local_symbol(2)
    tups = s2.symbol_tuple_list()
    for t in tups:
        n = t[1]
        d = t[2]
        s = t[3]
        o = t[4]
        if n == 0:
           if (s == 1) or (d in [3,5]):
               return False
        if n == 1:
            if (d in [1,7]) and (o not in [1,7]):
                return False
            if (d in [3,5]) and (o not in [3,5]):
                return False
        if (n == 2) and (s == 1):
            if (d in [1,7]) and (o not in [0,4,6]):
                return False
            if (d in [3,5]) and (o not in [2,4,6]):
                return False
        if (n > 2):
            if ((o - n) % 2 != 0):
                return False
            if (s == 0) and (o != 0):
                return False
    return True

def check_genus_symbol(genus_symbol):
    jordan = jordan_condition(genus_symbol)
    odd = oddity_condition(genus_symbol)
    ls = genus_symbol.local_symbols()
    det = genus_symbol.determinant()
    det_conds = [(sp.prime(), determinant_condition(sp, det)) for sp in ls]
    return det_conds, odd, jordan

# (Working on code for filling all the entries in lat_lattices and lat_genera)
