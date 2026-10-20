# This script takes input files from Noam Elkies' GP scripts and creates upload files for lat_lattices and lat_genera
# All functions in this script should be run from this folder (since various paths are relative); note that this will probably require doing sys.path.append(LOCATION_OF_LMFDB) first.

magma.attach("canonical_form.m")
from collections import defaultdict
from sage.databases.cremona import class_to_int, cremona_letter_code
from lmfdb.backend.encoding import copy_dumps
import re
opj = os.path.join
ope = os.path.exists
md_link_re = re.compile(r"\s*\[(.*)\]\((.*)\)\s*")

def parse_line(line):
    invs, gram = line.strip().split("][")
    res = {}
    res["gram"] = sage_eval("[" + gram.replace(";", ","))
    res["root_string"], res["mw_rank"], res["mw_torsion"], res["root_discriminant"], res["num_reducible_fibers"], res["theta_series"], res["quest"] = sage_eval(invs + "]")
    res["root_lattice"] = parse_root_string(res["root_string"])
    return res

def parse_root_string(s):
    pieces = []
    for part in s.split():
        if "^" not in part:
            n = 1
        else:
            part, n = part.split("^")
            n = ZZ(n)
        for i in range(n):
            pieces.append((part[0], ZZ(part[1:])))
    pieces.sort(reverse=True)
    return tuple(pieces)

def process_genus(genus_label):
    """
    Process a positive definite genus that has been enumerated by Noam's code
    """
    infile = opj("data", "nde_out", genus_label)
    ZZx = ZZ['x']
    with open(infile) as F:
        chunks = F.read().strip().split("\n\n")
        assert len(chunks) == 5
        lats = [parse_line(line) for line in chunks[0].split("\n")]
        genus_data = {}
        n = genus_data["class_number"] = len(lats)
        assert n-1 == chunks[1].count("\n") == chunks[2].count(";") == chunks[3].count("\n")
        for res, aut in zip(lats, chunks[1].split("\n")):
            res["aut_size"] = ZZ(aut)
        genus_data["adjacency_matrix"] = {"2": sage_eval(chunks[2].replace(";", ","))}
        polys = chunks[3].split("\n")
        assert all(f[0] == "[" and f[-1] == "]" for f in polys)
        polys = [f[1:-1].rsplit(" ", 1) for f in polys]
        polys = [[list(ZZx(f[0])), ZZ(f[1])] for f in polys]
        genus_data["adjacency_polynomials"] = {"2": polys}
        mass = chunks[4].lstrip("mass = ")
        if "/" in mass:
            mass = [ZZ(c) for c in mass.split("/")]
        else:
            mass = [ZZ(mass), 1]
        genus_data["mass"] = mass

    def sort_key(res):
        return (res["root_lattice"], prod(res["mw_torsion"]), tuple(res["mw_torsion"]), tuple(res["theta_series"]), res["quest"])

    by_skey = defaultdict(list)
    for lat in lats:
        by_skey[sort_key(lat)].append(lat)
    slats = []
    for key in sorted(by_skey, reverse=True):
        L = by_skey[key]
        if len(L) > 1:
            print("Computing canonical forms", key)
            for lat in L:
                M = magma.MatrixAlgebra(ZZ, ZZ(len(lat["gram"])).isqrt())
                A = M(lat["gram"])
                can = A.CanonicalForm()
                lat["canonical"] = [ZZ(can[i][j]) for i in range(1,ZZ(can.Nrows())+1) for j in range(1,ZZ(can.Ncols())+1)]
            L.sort(key=lambda lat: lat["canonical"], reverse=True)
        slats.extend(L)

    for i, lat in enumerate(slats):
        lat["label"] = f"{genus_label}.{i+1}"

    return genus_data, slats

def load_schema_hashes():
    fbase = opj("..", "schemas")
    old_hashes = {}
    if ope(opj(fbase, "hashes")):
        with open(opj(fbase, "hashes")) as F:
            for line in F:
                table, shash, dhash = line.strip().split()
                old_hashes[table] = (int(shash), int(dhash))
    return old_hashes

def save_schema_hashes(new_hashes):
    fbase = opj("..", "schemas")
    with open(opj(fbase, "hashes"), "w") as F:
        for table, (shash, dhash) in new_hashes.items():
            _ = F.write(f"{table} {shash} {dhash}\n")

def load_schemas(reset=False):
    """
    Load schema files from Markdown into a dictionary ready for being written
    """
    fbase = opj("..", "schemas")
    old_hashes = load_schema_hashes()
    schemas = defaultdict(list)
    descriptions = defaultdict(list)
    new_hashes = {}
    warned = False
    for fname in os.listdir(fbase):
        if fname.endswith(".md"):
            table = fname[:-3]
            with open(opj(fbase, fname)) as F:
                for line in F:
                    if line.count("|") == 4:
                        col, typ, desc = line.split("|")[1:-1]
                        col = col.strip()
                        typ = typ.strip()
                        # desc is not stripped above, so compare stripped here or
                        # the header row parses as a column literally named "Column".
                        if col == "Column" and typ == "Type" and desc.strip() == "Description":
                            continue
                        hline = set("-:")
                        if all(set(x.strip()).issubset(hline) for x in [col, typ, desc]):
                            continue
                        m = md_link_re.fullmatch(col)
                        if m:
                            col = m.group(1)
                        schemas[table].append((col, typ))
                        descriptions[table].append(desc)
            # tuple(), because schemas[table]/descriptions[table] are lists and
            # hash() rejects those.
            new_hashes[table] = (hash(tuple(schemas[table])), hash(tuple(descriptions[table])))
            if table in old_hashes:
                if old_hashes[table] != new_hashes[table]:
                    warned = True
                    if old_hashes[table][0] != new_hashes[table][0]:
                        wtype = "Columns/types"
                    else:
                        wtype = "Descriptions"
                    print(f"Warning!  {wtype} for {table} have changed, so stored data may not be valid.  Call load_schemas(reset=True) to suppress this message, and delete_stored_data('{table}') to delete all stored data for this table.")
    if reset or not old_hashes or (not warned and sorted(new_hashes) != sorted(old_hashes)):
        save_schema_hashes(new_hashes)
    return schemas, warned, descriptions

def write_location(table, label=None):
    if label is None:
        return opj("data", table)
    return opj("data", table, label)

def collate_location(table):
    return opj("data", f"{table}.txt")

def delete_stored_data(table=None):
    if table is None:
        for table in load_schema_hashes():
            delete_stored_data(table)
    else:
        os.unlink(write_location(table))

def ensure_folders(schemas, genus_label=None, overwrite=False):
    """
    Make sure that all data folders exist
    """
    for table in schemas:
        os.makedirs(write_location(table))
        if genus_label is not None:
            fname = write_location(table, genus_label)
            if not overwrite and ope(fname):
                raise ValueError(f"Would overwrite {fname}; call with overwrite=True to proceed")

def write_header(schema):
    return "|".join(col for (col, typ) in schema) + "\n" + "|".join(typ for (col, typ) in schema) + "\n\n"

def write_line(record, schema):
    return "|".join(copy_dumps(record.get(col), typ) for (col, typ) in schema) + "\n"

def write_upload_files(genus_label, overwrite=False, reset_schema_hashes=False):
    posdef_genus, lats = process_genus(genus_label)
    schemas, warned, _ = load_schemas(reset_schema_hashes)
    if warned and not overwrite:
        raise ValueError("Either reset schema hashes or use overwrite=True to proceed")
    ensure_folders(schemas, genus_label, overwrite)

    # From this data, we need to produce
    # 1 positive-definite entry in the lat_genera table, and 1 for the direct sum with the hyperbolic plane
    # class_number entries in the lat_lattices table for the positive definite latties, and 1 for the indefinite
    # 1 entry in the k3_families table
    # class_number entries in the k3_elliptic table
    # class_number entries in the k3_family_models table (the Weierstrass models of the elliptic surfaces)

    # Needed (easier):
    #  The genus label for L+U (and more if this is not unique in its genus)
    #  An explicit embedding of L+U+T into II_{3,19}
    #  Label of the transcendental lattice


    for table, data in [("lat_genera", posdef_genus),
                        ("lat_genera", indef_genus),
                        ("lat_lattices", lats),
                        ("lat_lattices", indef_lat),
                        ("k3_families", k3fam),
                        ("k3_elliptic", k3ells),
                        ("k3_family_models", ell_models)]:
        schema = schemas[table]
        with open(write_location(table, genus_label), "a") as F:
            if isinstance(data, list):
                for X in data:
                    _ = F.write(write_line(X, schema))
            else:
                _ = F.write(write_line(data, schema))

def collate_data_files():
    schemas, warned, _ = load_schemas()
    if warned:
        raise ValueError("Reset schema hashes to proceed")
    ensure_folders(schemas)
    for table, schema in schemas.items():
        with open(collate_location(table), "w") as Fout:
            _ = Fout.write(write_header(schema))
            for fname in os.listdir(write_location(table)):
                with open(write_location(table, fname)) as F:
                    for line in F:
                        _ = Fout.write(line)

def create_tables():
    schemas, warned, descriptions = load_schemas()
    if warned:
        raise ValueError("Reset schema hashes to proceed")
    for table, schema in schemas.items():
        search_columns = defaultdict(list)
        for col, typ in schema:
            search_columns[typ].append(col)
        search_order = [col for (col, typ) in schema]
        col_description = dict(zip(search_order, descriptions[table]))
        # TODO: Is this always right?
        label_col = "label" if "label" in col_description else None
        # TODO: Ugh, this should be added to schema
        table_description = table
        # TODO: This should be added to the schema
        sort = None
        db.create_table(table,
                        seach_columns,
                        label_col,
                        table_description=table_description,
                        col_description=col_description,
                        sort=sort,
                        search_order=search_order)
