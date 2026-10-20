"""
Check schemas/*.md against the .format files the pipeline actually writes.
Run: python3 tests/test_schema_consistency.py   (from lattices/; no Sage needed)

Why this exists: schemas/*.md is not just documentation.  nde_to_upload.sage's
load_schemas() parses those markdown tables into the (column, type) list that
write_header()/write_line() use, and write_line() does record.get(col) -- so a
column name in the markdown is a lookup key into the record dict.  Get it wrong
and the upload silently writes nulls instead of failing.  Five such bugs were
live when this check was written, none visible to the eye:

  * "| deep_hole_count | " had only 2 pipes, so the parser skipped the row
    entirely and the column vanished from the upload.
  * theta_prec appeared twice (the second meant dual_theta_prec).  Schemas
    accumulate in a list, not a dict, so the upload got two identical
    theta_prec columns and no dual_theta_prec at all.
  * covering_norm was one row, but the pipeline emits covering_norm_num and
    covering_norm_den.
  * "[canonical_gram]" had brackets with no link URL, so md_link_re (which
    requires "[...](...)") did not strip them and the column parsed as
    "[canonical_gram]", brackets included.
  * the header row was skipped by comparing an UNSTRIPPED description against
    "Description", which never matched, so every table gained a leading column
    literally named "Column".

That last one is why this file drives the real load_schemas() -- extracted from
nde_to_upload.sage and exec'd -- instead of reimplementing the parse.  An earlier
version of this check mirrored the parser by hand, and the copy was *more*
correct than the original, which hid the bug.  Only the three schema functions
are extracted, so this needs no Sage and no lmfdb import.
"""
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
LATTICES = os.path.join(HERE, "..")

# Which .format files together make up each table's columns.
TABLE_FORMATS = {
    "lat_genera": ["genera_basic", "genera_advanced"],
    "lat_lattices": ["lat_basic", "lat_advanced", "lat_decomp", "lat_hash"],
}

# Columns intentionally declared in the schema before anything computes them.
# Removing one from this set should mean the pipeline now populates it.
NOT_YET_COMPUTED = {
    "lat_lattices": {"even_complement"},
}


def real_load_schemas():
    """Run nde_to_upload.sage's own load_schemas(), without importing the file.

    Importing it outright would pull in sage.databases.cremona and
    lmfdb.backend.encoding; the three schema functions only need os/re/
    defaultdict, so extract just those.
    """
    src = open(os.path.join(LATTICES, "nde_to_upload.sage")).read()

    def grab(name):
        start = src.index("def %s(" % name)
        end = src.index("\ndef ", start + 1)
        return src[start:end]

    ns = {
        "os": os,
        "re": re,
        "defaultdict": defaultdict,
        "opj": os.path.join,
        "ope": os.path.exists,
        "md_link_re": re.compile(r"\s*\[(.*)\]\((.*)\)\s*"),
    }
    for fn in ("load_schema_hashes", "save_schema_hashes", "load_schemas"):
        exec(grab(fn), ns)
    # Never touch schemas/hashes from a test run.
    ns["save_schema_hashes"] = lambda new_hashes: None

    cwd = os.getcwd()
    try:
        os.chdir(LATTICES)   # load_schemas() resolves schemas/ as ../schemas
        schemas, _warned, _descriptions = ns["load_schemas"]()
    finally:
        os.chdir(cwd)
    return schemas


def suspicious_rows(table):
    """Rows that look like table rows but that load_schemas() would drop."""
    path = os.path.join(LATTICES, "..", "schemas", table + ".md")
    out = []
    for lineno, line in enumerate(open(path), 1):
        if line.lstrip().startswith("|") and line.count("|") != 4:
            out.append((lineno, line.rstrip()))
    return out


def format_columns(name):
    with open(os.path.join(LATTICES, name + ".format")) as f:
        header = f.read().split("\n")[0]
    return [c.strip() for c in header.split("|")]


def main():
    schemas = real_load_schemas()
    failures = []

    for table, fmt_names in sorted(TABLE_FORMATS.items()):
        names = [col for col, _typ in schemas[table]]

        for lineno, text in suspicious_rows(table):
            failures.append(
                "%s.md:%d: row would be silently dropped by load_schemas() "
                "(needs exactly 4 '|'): %r" % (table, lineno, text)
            )

        for name in sorted({n for n in names if names.count(n) > 1}):
            failures.append(
                "%s.md: column %r parsed %d times; schemas are a list, so this "
                "duplicates a column in the upload" % (table, name, names.count(name))
            )

        for col, typ in schemas[table]:
            if col != col.strip("[]"):
                failures.append(
                    "%s.md: column %r still has brackets -- only a full "
                    "'[text](url)' link gets stripped" % (table, col)
                )
            if not typ:
                failures.append("%s.md: column %r has no type" % (table, col))
            if col in ("Column", "Type", "Description"):
                failures.append(
                    "%s.md: parsed a column named %r -- the markdown header row "
                    "is being treated as data" % (table, col)
                )

        produced = set()
        for fmt_name in fmt_names:
            produced |= set(format_columns(fmt_name))

        allowed = NOT_YET_COMPUTED.get(table, set())
        for name in sorted(set(names) - produced - allowed):
            failures.append(
                "%s.md: column %r is declared but no .format file (%s) produces "
                "it -- the upload would write nulls for it. If that is intended, "
                "add it to NOT_YET_COMPUTED." % (table, name, ", ".join(fmt_names))
            )
        for name in sorted(produced - set(names)):
            failures.append(
                "%s: column %r is written by the pipeline but is absent from "
                "%s.md, so it would never be uploaded" % (table, name, table)
            )
        for name in sorted(allowed & produced):
            failures.append(
                "%s: column %r is in NOT_YET_COMPUTED but the pipeline now "
                "produces it -- remove it from that set" % (table, name)
            )

        print("  %s: %d columns, %d produced" % (table, len(names), len(produced)))

    print()
    if failures:
        print("FAIL: %d schema/pipeline inconsistencies" % len(failures))
        for f in failures:
            print("  - %s" % f)
        sys.exit(1)
    print("PASS: schemas and .format files agree")


if __name__ == "__main__":
    main()
