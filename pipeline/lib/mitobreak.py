#!/usr/bin/env python3
"""MitoBreak known-breakpoint database: load + match a called SV to the closest
previously-reported mtDNA breakpoint.

Mirrors MitoHPC's svMitoBreak annotation so the cohort report can count KNOWN
(catalogued in MitoBreak) vs NOVEL calls per caller — a credibility/specificity
proxy a reviewer can read at a glance.

DB: assets/mitobreak.tsv.gz, columns
  svtype  bp5  bp3  length  location  origin_impact  disease  refs  mitobreak_id
1-based rCRS (NC_012920.1). Vendored verbatim from MitoHPC RefSeq/mitobreak.tsv.gz
(1369 deletions + 44 duplications).
"""
from __future__ import annotations

import csv
import gzip
import os

# Per-breakpoint tolerance (bp). MitoHPC's svMitoBreak default is 20; we keep it —
# it absorbs the ~13 bp direct-repeat spread of the common deletion and the +/-1
# breakpoint-convention jitter between callers, while keeping spurious hits on the
# dense catalogue rare. Applied identically to every caller, so the known/novel
# split is a fair cross-caller comparison.
MITOBREAK_TOL = 20

# our normalised sv_type -> MitoBreak svtype (the DB only catalogues DEL and DUP)
_TYPE = {"deletion": "DEL", "duplication": "DUP"}


def _int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def default_db():
    """Locate the vendored DB: the image copy first, then the repo copy."""
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ("/opt/assets/mitobreak.tsv.gz",
              os.path.join(here, "..", "..", "assets", "mitobreak.tsv.gz")):
        if os.path.isfile(p):
            return p
    return None


def load(path=None):
    """-> {svtype: [(bp5, bp3, length, mid, disease, location), ...]} (empty if absent)."""
    path = path or default_db()
    db = {}
    if not path or not os.path.isfile(path):
        return db
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            b5, b3 = _int(row.get("bp5")), _int(row.get("bp3"))
            t = (row.get("svtype") or "").strip()
            if b5 is None or b3 is None or not t:
                continue
            db.setdefault(t, []).append(
                (b5, b3, _int(row.get("length")),
                 (row.get("mitobreak_id") or "").strip(),
                 (row.get("disease") or "").strip(),
                 (row.get("location") or "").strip()))
    return db


def match(db, sv_type, bp5, bp3, tol=MITOBREAK_TOL):
    """Closest MitoBreak entry whose 5' and 3' breakpoints are BOTH within `tol` of
    a same-type call, or None. Returns (mid, disease, location, summed_dist)."""
    if bp5 is None or bp3 is None:
        return None
    t = _TYPE.get(sv_type)
    rows = db.get(t) if t else None
    if not rows:
        return None
    # our DEL bp3 is the last-DELETED base; MitoBreak's bp3 is the first-RETAINED
    # base after the event (exclusive). +1 reconciles the conventions (MitoHPC does
    # the same); the difference is < tol anyway, so it never changes a match.
    q3 = bp3 + 1 if t == "DEL" else bp3
    best = None
    for (m5, m3, _mlen, mid, dis, loc) in rows:
        d5, d3 = abs(bp5 - m5), abs(q3 - m3)
        if d5 <= tol and d3 <= tol:
            d = d5 + d3
            if best is None or d < best[0]:
                best = (d, mid, dis, loc)
    return None if best is None else (best[1], best[2], best[3], best[0])
