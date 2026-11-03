#!/usr/bin/env python3
"""Assemble the pipeline's flat files into LMFDB upload tables, checking every value.

    python3 make_upload.py [--out upload] [--limit N]

Writes <out>/lat_genera.txt and <out>/lat_lattices.txt in psycodict's reload format
(column names, column types, a blank line, then one '|'-separated row per record),
with columns in schemas/*.md order, and copies lat_inclusions.tbl if present.
Prints, per table, the null count of every column and every value that does not fit
its column type; exits 1 if there are any (nothing that would fail or silently
change in Postgres goes unreported).

Sources (run from lattices/ after connect):
  lat_genera    genera_advanced/<r>/<nplus>/<label>  (genera_basic + genera_advanced
                fields); a listed genus without one (never filled) gets its
                genera_basic fields and nulls
  lat_lattices  lattice_basic_data + lattice_advanced_data (per label),
                lattice_decomp_data/<rank>.<nplus> (tensor stage) and the per-genus
                lattice_hashes files, joined on the label

The schema is the one nde_to_upload.sage's load_schemas() parses (run here through
tests/test_schema_consistency.py, so there is a single parser).  This replaces no
upload script: nde_to_upload.sage handles a different input (K3 families from
data/nde_out/) and never read these files.
"""
import argparse
import json
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests"))
from test_schema_consistency import real_load_schemas  # noqa: E402

NULL = "\\N"
INT_RANGE = {"smallint": 2**15, "integer": 2**31, "bigint": 2**63}


def read_format(name):
    return Path(f"{name}.format").read_text().strip().split("|")


def split_array(s):
    """Top-level elements of a Postgres array literal {a,b,{c,d}}."""
    inner = s[1:-1]
    out, depth, cur = [], 0, []
    for ch in inner:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if cur or out:
        out.append("".join(cur))
    return [x.strip() for x in out]


def check_scalar(v, typ):
    """None if the text v is a valid typ value (as Postgres COPY reads it), else why."""
    if typ in INT_RANGE:
        if not re.fullmatch(r"-?\d+", v):
            return "not an integer"
        if not -INT_RANGE[typ] <= int(v) < INT_RANGE[typ]:
            return f"out of {typ} range"
        return None
    if typ == "numeric":
        try:
            Decimal(v)
        except InvalidOperation:
            return "not a numeric literal" + (" (a fraction)" if "/" in v else "")
        return None
    if typ in ("float4", "float8", "real", "double precision"):
        try:
            float(v)
        except ValueError:
            return "not a float"
        return None
    if typ == "boolean":
        return None if v in ("t", "f") else "not a boolean"
    if typ == "text":
        return None
    if typ == "jsonb":
        try:
            json.loads(v)
        except ValueError:
            return "not JSON"
        return None
    return f"unknown type {typ}"


def check_value(v, typ):
    if v == NULL:
        return None
    if typ.endswith("[]"):
        if not (v.startswith("{") and v.endswith("}")):
            return "not an array"
        for x in split_array(v):
            if x.startswith("{"):
                why = check_value(x, typ)          # a nested array: already says "element"
                if why:
                    return why
            elif x != "NULL":
                why = check_scalar(x, typ[:-2])
                if why:
                    return "element " + why
        return None
    return check_scalar(v, typ)


def normalize(v, typ):
    """Pipeline text -> COPY text for this type."""
    if v in ("None", ""):
        return NULL
    if typ == "boolean":
        return {"T": "t", "F": "f", "true": "t", "false": "f"}.get(v, v)
    if typ == "text":
        # COPY text format: escape backslashes and the delimiter; \N stays a null
        return v if v == NULL else v.replace("\\", "\\\\").replace("|", "\\|").replace("\n", "\\n")
    return v


class Table:
    def __init__(self, name, schema):
        self.name, self.schema = name, schema
        self.nulls = Counter()
        self.bad = defaultdict(list)          # column -> [(label, value, why)]
        self.rows = 0

    def row(self, label, rec):
        out = []
        for col, typ in self.schema:
            v = normalize(rec.get(col, NULL), typ)
            if v == NULL:
                self.nulls[col] += 1
            why = check_value(v, typ)
            if why:
                self.bad[col].append((label, v[:80], why))
            out.append(v)
        self.rows += 1
        return "|".join(out) + "\n"

    def header(self):
        return "|".join(c for c, _ in self.schema) + "\n" + "|".join(t for _, t in self.schema) + "\n\n"

    def report(self, produced):
        print(f"\n{self.name}: {self.rows} rows")
        cols = {c for c, _ in self.schema}
        extra = sorted(set(produced) - cols)
        if extra:
            print(f"  produced but not in the schema (dropped): {', '.join(extra)}")
        never = [c for c, _ in self.schema if c not in produced]
        if never:
            print(f"  in the schema but produced by no stage (all null): {', '.join(never)}")
        for col, typ in self.schema:
            k = self.nulls[col]
            if k and col not in never:
                print(f"  {col:32s} {typ:10s} null {k:8d} ({k / max(self.rows, 1):6.1%})")
        nbad = sum(len(v) for v in self.bad.values())
        for col, items in self.bad.items():
            print(f"  INVALID {col}: {len(items)} values, e.g. " +
                  "; ".join(f"{l}: {v!r} ({w})" for l, v, w in items[:3]))
        return nbad


def label_files(folder):
    for dp, _, fns in os.walk(folder):
        for fn in fns:
            if not fn.endswith((".tmp", ".done")) and ".tmp." not in fn:
                yield fn, os.path.join(dp, fn)


def first_line_fields(path):
    with open(path) as F:
        return F.readline().rstrip("\n").split("|")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="upload")
    ap.add_argument("--limit", type=int, default=0, help="stop after this many genera (for a quick check)")
    a = ap.parse_args()
    schemas = real_load_schemas()
    out = Path(a.out)
    out.mkdir(exist_ok=True)
    gb, ga = read_format("genera_basic"), read_format("genera_advanced")
    lb, la, ld, lh = (read_format(n) for n in ("lat_basic", "lat_advanced", "lat_decomp", "lat_hash"))
    nbad = 0

    # lat_genera
    T = Table("lat_genera", schemas["lat_genera"])
    produced = set(gb) | set(ga)
    filled = 0
    with open(out / "lat_genera.txt", "w") as F:
        F.write(T.header())
        for n, (label, bpath) in enumerate(label_files("genera_basic")):
            if a.limit and n >= a.limit:
                break
            apath = os.path.join("genera_advanced", *bpath.split(os.sep)[-3:])
            if os.path.exists(apath):
                fields = first_line_fields(apath)
                rec = dict(zip(gb + ga, fields)) if len(fields) == len(gb + ga) else None
                filled += 1
            else:
                rec = None
            if rec is None:
                fields = first_line_fields(bpath)
                rec = dict(zip(gb, fields))
            F.write(T.row(label, rec))
    nbad += T.report(produced)
    print(f"  {filled} of {T.rows} genera filled")

    # lat_lattices: decompositions and hashes are keyed by label in bulk files
    decomp = {}
    if os.path.isdir("lattice_decomp_data"):
        for _, p in label_files("lattice_decomp_data"):
            for line in open(p):
                f = line.rstrip("\n").split("|")
                if len(f) == len(ld):
                    decomp[f[0]] = dict(zip(ld, f))
    hashes = {}
    for _, p in label_files("lattice_hashes"):
        for line in open(p):
            f = line.rstrip("\n").split("|")
            if len(f) == len(lh):
                hashes[f[0]] = dict(zip(lh, f))
    T = Table("lat_lattices", schemas["lat_lattices"])
    produced = set(lb) | set(la) | set(ld) | set(lh)
    no_adv = 0
    with open(out / "lat_lattices.txt", "w") as F:
        F.write(T.header())
        for label, bpath in label_files("lattice_basic_data"):
            rec = dict(zip(lb, first_line_fields(bpath)))
            apath = os.path.join("lattice_advanced_data", *bpath.split(os.sep)[-3:])
            if os.path.exists(apath):
                f = first_line_fields(apath)
                if len(f) == len(la):
                    rec.update(dict(zip(la, f)))
            else:
                no_adv += 1
            rec.update(decomp.get(label, {}))
            rec.update({k: v for k, v in hashes.get(label, {}).items() if k != "gram"})
            F.write(T.row(label, rec))
    nbad += T.report(produced)
    if no_adv:
        print(f"  {no_adv} lattices have no advanced record (connect did not run on them)")
    if os.path.exists("lat_inclusions.tbl"):
        shutil.copyfile("lat_inclusions.tbl", out / "lat_inclusions.txt")
        print("\nlat_inclusions: copied lat_inclusions.tbl (checked by oscar/write_table.jl)")
    print(f"\n{nbad} invalid values" + ("" if nbad == 0 else " -- fix before uploading"))
    return 1 if nbad else 0


if __name__ == "__main__":
    sys.exit(main())
