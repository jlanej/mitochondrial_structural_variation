"""Parsers that normalise each caller's SV output into a common schema.

Every parser returns a list of ``dict`` records with these keys:

    sample   str    sample identifier
    caller   str    one of eklipse|mitosalt|splicebreak2|mitomut|mitoseek
    sv_type  str    deletion|duplication|breakpoint|sv
    bp5      int    5' breakpoint (1-based, rCRS), or None
    bp3      int    3' breakpoint (1-based, rCRS), or None
    svlen    int    event size in bp, or None
    support  float  number of supporting reads, or None
    het      float  heteroplasmy / allele fraction in [0,1], or None
    extra    str    free-text caller-specific detail

Coordinates are reported on the canonical rCRS (NC_012920.1 / "chrM", 16569 bp).
Heteroplasmy is always normalised to a 0-1 fraction.

These functions are deliberately defensive: a malformed row is skipped, never
fatal, so one ragged line never sinks a whole cohort merge.
"""
from __future__ import annotations

import csv
import os
import re
from statistics import median

CALLERS = ["eklipse", "mitosalt", "splicebreak2", "mitomut", "mitoseek"]

# The classic ~4977 bp "common deletion" (m.8470_13447del). Callers place the
# breakpoints anywhere inside the flanking 13 bp direct repeat, so we match with
# tolerance.
COMMON_DEL_BP5 = 8470
COMMON_DEL_BP3 = 13447
COMMON_DEL_TOL = 80


def _rec(sample, caller, sv_type="deletion", bp5=None, bp3=None,
         svlen=None, support=None, het=None, extra=""):
    return {
        "sample": sample, "caller": caller, "sv_type": sv_type,
        "bp5": bp5, "bp3": bp3, "svlen": svlen,
        "support": support, "het": het, "extra": extra,
    }


def _to_int(x):
    try:
        return int(round(float(str(x).strip())))
    except (ValueError, TypeError):
        return None


def _to_float(x):
    try:
        return float(str(x).strip().replace(",", "."))
    except (ValueError, TypeError):
        return None


def is_common_deletion(bp5, bp3):
    """True if (bp5, bp3) matches the ~4977 bp common deletion within tolerance."""
    if bp5 is None or bp3 is None:
        return False
    return (abs(bp5 - COMMON_DEL_BP5) <= COMMON_DEL_TOL
            and abs(bp3 - COMMON_DEL_BP3) <= COMMON_DEL_TOL)


# --------------------------------------------------------------------------- #
# eKLIPse — eKLIPse_deletions.csv  (";"-delimited, quoted, comma decimals)
# --------------------------------------------------------------------------- #
def parse_eklipse(path, sample=None):
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        for i, raw in enumerate(fh):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = [f.strip().strip('"') for f in line.split(";")]
            if i == 0 or fields[0].lower().startswith("title"):
                continue  # header
            if len(fields) < 3:
                continue
            title = fields[0]
            bp5 = _to_int(fields[1])
            bp3 = _to_int(fields[2])
            freq = _to_float(fields[3]) if len(fields) > 3 else None
            d5 = _to_float(fields[8]) if len(fields) > 8 else None
            d3 = _to_float(fields[9]) if len(fields) > 9 else None
            rep = fields[10] if len(fields) > 10 else ""
            svlen = (bp3 - bp5) if (bp5 is not None and bp3 is not None) else None
            het = (freq / 100.0) if freq is not None else None
            support = min(x for x in (d5, d3) if x is not None) if (d5 or d3) else None
            out.append(_rec(sample or title, "eklipse", "deletion",
                            bp5, bp3, svlen, support, het,
                            extra=("repeat=" + rep) if rep else ""))
    return out


# --------------------------------------------------------------------------- #
# MitoSAlt — indel/<tag>.tsv  (tab, header row)
# --------------------------------------------------------------------------- #
def parse_mitosalt(path, sample=None):
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = [r for r in reader if r]
    if not rows:
        return out
    header = [h.strip() for h in rows[0]]
    idx = {h: k for k, h in enumerate(header)}

    def col(row, *names):
        for n in names:
            if n in idx and idx[n] < len(row):
                return row[idx[n]]
        return None

    for row in rows[1:]:
        if not row or row[0].startswith("#"):
            continue
        ev = (col(row, "final.event", "event") or "deletion").strip().lower()
        sv_type = "duplication" if "dup" in ev else "deletion"
        bp5 = _to_int(col(row, "final.start", "del.start", "start"))
        bp3 = _to_int(col(row, "final.end", "del.end", "end"))
        svlen = _to_int(col(row, "final.size", "delsize", "size"))
        het = _to_float(col(row, "heteroplasmy"))
        support = _to_float(col(row, "alt.reads", "alt"))
        smp = sample or (col(row, "sample") or "").strip() or "unknown"
        # delplot.R names the direct-repeat flank column "seq" (older builds:
        # "direct.repeat"); final.event values are "del"/"dup".
        rep = (col(row, "seq", "direct.repeat") or "").strip()
        out.append(_rec(smp, "mitosalt", sv_type, bp5, bp3, svlen, support, het,
                        extra=("repeat=" + rep) if rep else ""))
    return out


# --------------------------------------------------------------------------- #
# Splice-Break2 — *_LargeMTDeletions_*.txt  (whitespace-aligned, header row)
# --------------------------------------------------------------------------- #
def parse_splicebreak2(path, sample=None):
    out = []
    if not path or not os.path.isfile(path):
        return out
    header = None
    idx = {}
    with open(path) as fh:
        for line in fh:
            toks = line.split()
            if not toks:
                continue
            if header is None:
                if toks[0] == "Sample_ID":
                    header = toks
                    idx = {h: k for k, h in enumerate(header)}
                continue
            # Data row.
            def col(name):
                k = idx.get(name)
                return toks[k] if (k is not None and k < len(toks)) else None
            smp = sample or col("Sample_ID") or "unknown"
            bp5 = _to_int(col("5'_Break"))
            bp3 = _to_int(col("3'_Break"))
            svlen = _to_int(col("Deletion_Size_bp"))
            support = _to_float(col("Deletion_Reads"))
            pct = _to_float(col("Deletion_Read_%"))
            het = (pct / 100.0) if pct is not None else None
            ann = col("Annotation") or ""
            out.append(_rec(smp, "splicebreak2", "deletion", bp5, bp3, svlen,
                            support, het, extra=("annotation=" + ann) if ann else ""))
    return out


# --------------------------------------------------------------------------- #
# MitoMut — <bam>_results.txt  (tab, header row)
# --------------------------------------------------------------------------- #
def parse_mitomut(path, sample=None):
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = [r for r in reader if r]
    for row in rows:
        if not row or row[0].strip().lower().startswith("total"):
            continue  # header
        if len(row) < 6:
            continue
        total, s1, s2, start, end, het = row[:6]
        bp5 = _to_int(start)
        bp3 = _to_int(end)
        svlen = (abs(bp3 - bp5) if (bp5 is not None and bp3 is not None) else None)
        out.append(_rec(sample or "unknown", "mitomut", "deletion",
                        bp5, bp3, svlen, _to_float(total), _to_float(het),
                        extra="s1=%s;s2=%s" % (s1.strip(), s2.strip())))
    return out


# --------------------------------------------------------------------------- #
# MitoSeek — discordant mates table + large-deletion SAM
# --------------------------------------------------------------------------- #
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _ref_len_from_cigar(cigar):
    if not cigar or cigar == "*":
        return 0
    return sum(int(n) for n, op in _CIGAR_RE.findall(cigar) if op in "MDN=X")


def parse_mitoseek_discordant(path, sample=None):
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            mitopos = _to_int(f[1])
            support = _to_float(f[4])
            out.append(_rec(sample or "unknown", "mitoseek", "breakpoint",
                            bp5=mitopos, bp3=None, svlen=None, support=support,
                            het=None, extra="mate=%s:%s" % (f[2], f[3])))
    return out


def parse_mitoseek_large_deletion(path, sample=None, min_tlen=500, min_support=3,
                                  bin_size=25):
    """Aggregate large-deletion supporting reads in the MitoSeek SAM into calls.

    MitoSeek emits the raw read pairs whose template length exceeds -strf rather
    than a summarised call set, so we cluster them by approximate breakpoints.
    """
    out = []
    if not path or not os.path.isfile(path):
        return out
    clusters = {}  # (bp5_bin, bp3_bin) -> list of (bp5, bp3, svlen)
    with open(path) as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            flag = _to_int(f[1]) or 0
            # primary, mapped, leftmost mate of the pair (positive TLEN)
            if flag & 0x100 or flag & 0x800 or flag & 0x4:
                continue
            tlen = _to_int(f[8])
            if tlen is None or tlen < min_tlen:
                continue
            pos = _to_int(f[3])
            pnext = _to_int(f[7])
            if pos is None or pnext is None:
                continue
            reflen = _ref_len_from_cigar(f[5])
            bp5 = pos + reflen          # end of left mate ~ 5' breakpoint
            bp3 = pnext                 # start of right mate ~ 3' breakpoint
            key = (bp5 // bin_size, bp3 // bin_size)
            clusters.setdefault(key, []).append((bp5, bp3, tlen))
    for key, members in clusters.items():
        if len(members) < min_support:
            continue
        bp5 = int(median([m[0] for m in members]))
        bp3 = int(median([m[1] for m in members]))
        svlen = bp3 - bp5
        out.append(_rec(sample or "unknown", "mitoseek", "deletion",
                        bp5, bp3, svlen, support=float(len(members)), het=None,
                        extra="span_reads=%d" % len(members)))
    return out


# --------------------------------------------------------------------------- #
# Sample-directory dispatch
# --------------------------------------------------------------------------- #
def _first(globs):
    import glob
    for g in globs:
        hits = sorted(glob.glob(g))
        if hits:
            return hits[0]
    return None


def parse_sample_dir(sample_dir, sample=None):
    """Parse every caller's output found under one sample directory."""
    sample = sample or os.path.basename(os.path.normpath(sample_dir))
    records = []

    ek = _first([os.path.join(sample_dir, "eklipse", "eKLIPse_deletions.csv")])
    records += parse_eklipse(ek, sample)

    ms = _first([os.path.join(sample_dir, "mitosalt", "*.mitosalt.tsv"),
                 os.path.join(sample_dir, "mitosalt", "work", "indel", "*.tsv")])
    records += parse_mitosalt(ms, sample)

    sb = _first([os.path.join(sample_dir, "splicebreak2",
                              "*_LargeMTDeletions_WGS-only_NoPositionFilter.txt")])
    records += parse_splicebreak2(sb, sample)

    mm = _first([os.path.join(sample_dir, "mitomut", "mitomut_results.txt"),
                 os.path.join(sample_dir, "mitomut", "*_results.txt")])
    records += parse_mitomut(mm, sample)

    disc = _first([os.path.join(sample_dir, "mitoseek", "mitoseek_discordant_mates.txt"),
                   os.path.join(sample_dir, "mitoseek", "*", "mito1_structure_discordant_mates.txt")])
    records += parse_mitoseek_discordant(disc, sample)
    dsam = _first([os.path.join(sample_dir, "mitoseek", "mitoseek_large_deletion.sam"),
                   os.path.join(sample_dir, "mitoseek", "*", "mito1_structure_large_deletion.sam")])
    records += parse_mitoseek_large_deletion(dsam, sample)

    return records
