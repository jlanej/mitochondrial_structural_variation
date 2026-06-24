#!/usr/bin/env python3
"""Consolidate per-sample caller outputs into a cohort SV summary.

Scans an output root that contains one sub-directory per sample, each holding
per-caller sub-directories (eklipse/ mitosalt/ splicebreak2/ mitomut/ mitoseek/),
parses every caller's output via ``parsers.py``, and writes:

  cohort_sv_calls.tsv        long table: every normalised call from every caller
  cohort_common_deletion.tsv per sample x caller: was the ~4977 common deletion
                             detected, with its coordinates
  cohort_caller_matrix.tsv   sample x caller matrix of deletion-call counts
  cohort_summary.txt         human-readable digest

Usage:
  postprocess.py --root <output-root> [--out-dir <dir>]
  postprocess.py --sample-dir <dir> [--sample NAME] [--out-dir <dir>]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import parsers  # noqa: E402  (path set above)

CALLERS = parsers.CALLERS
CALL_COLUMNS = ["sample", "caller", "sv_type", "bp5", "bp3", "svlen",
                "support", "het", "common_deletion", "extra"]


def _looks_like_sample_dir(path):
    if not os.path.isdir(path):
        return False
    if os.path.isfile(os.path.join(path, "status.tsv")):
        return True
    return any(os.path.isdir(os.path.join(path, c)) for c in CALLERS)


def discover_sample_dirs(root):
    out = []
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if _looks_like_sample_dir(p):
            out.append(p)
    return out


def _fmt(v):
    return "" if v is None else (str(v) if not isinstance(v, float) else ("%.4g" % v))


def write_calls(records, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(CALL_COLUMNS)
        for r in records:
            common = 1 if parsers.is_common_deletion(r["bp5"], r["bp3"]) else 0
            w.writerow([r["sample"], r["caller"], r["sv_type"],
                        _fmt(r["bp5"]), _fmt(r["bp3"]), _fmt(r["svlen"]),
                        _fmt(r["support"]), _fmt(r["het"]), common, r["extra"]])


def write_common_deletion(records, samples, path):
    """One row per (sample, caller): best common-deletion match if any."""
    best = {}  # (sample, caller) -> record
    for r in records:
        if not parsers.is_common_deletion(r["bp5"], r["bp3"]):
            continue
        key = (r["sample"], r["caller"])
        prev = best.get(key)
        if prev is None or (r["support"] or 0) > (prev["support"] or 0):
            best[key] = r
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "caller", "common_deletion_detected",
                    "bp5", "bp3", "svlen", "support", "het"])
        for s in samples:
            for c in CALLERS:
                r = best.get((s, c))
                if r:
                    w.writerow([s, c, 1, _fmt(r["bp5"]), _fmt(r["bp3"]),
                                _fmt(r["svlen"]), _fmt(r["support"]), _fmt(r["het"])])
                else:
                    w.writerow([s, c, 0, "", "", "", "", ""])


def write_matrix(records, samples, path):
    """sample x caller matrix of deletion-call counts (+ common-del caller count)."""
    counts = {(s, c): 0 for s in samples for c in CALLERS}
    common = {s: set() for s in samples}
    for r in records:
        if r["sv_type"] == "deletion":
            counts[(r["sample"], r["caller"])] = counts.get((r["sample"], r["caller"]), 0) + 1
        if parsers.is_common_deletion(r["bp5"], r["bp3"]):
            common[r["sample"]].add(r["caller"])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + CALLERS + ["total_deletion_calls",
                                           "n_callers_common_deletion"])
        for s in samples:
            row = [counts.get((s, c), 0) for c in CALLERS]
            w.writerow([s] + row + [sum(row), len(common[s])])


def write_summary(records, samples, path):
    by_caller = {c: 0 for c in CALLERS}
    for r in records:
        if r["sv_type"] == "deletion":
            by_caller[r["caller"]] = by_caller.get(r["caller"], 0) + 1
    with open(path, "w") as fh:
        fh.write("Mitochondrial SV cohort summary\n")
        fh.write("=" * 40 + "\n")
        fh.write("samples: %d\n" % len(samples))
        fh.write("total calls (all types): %d\n" % len(records))
        fh.write("\ndeletion calls per caller:\n")
        for c in CALLERS:
            fh.write("  %-14s %d\n" % (c, by_caller.get(c, 0)))
        fh.write("\ncommon (~4977 bp) deletion detection per sample:\n")
        common = {s: set() for s in samples}
        for r in records:
            if parsers.is_common_deletion(r["bp5"], r["bp3"]):
                common[r["sample"]].add(r["caller"])
        for s in samples:
            cs = sorted(common[s])
            fh.write("  %-24s %s\n" % (s, ", ".join(cs) if cs else "(none)"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--root", help="output root containing per-sample sub-dirs")
    g.add_argument("--sample-dir", help="a single per-sample directory")
    ap.add_argument("--sample", help="sample name (with --sample-dir)")
    ap.add_argument("--out-dir", help="where to write cohort_*.tsv (default: root)")
    args = ap.parse_args(argv)

    if args.root:
        root = os.path.abspath(args.root)
        sample_dirs = discover_sample_dirs(root)
        out_dir = args.out_dir or root
    else:
        sd = os.path.abspath(args.sample_dir)
        sample_dirs = [sd]
        out_dir = args.out_dir or os.path.dirname(sd)

    os.makedirs(out_dir, exist_ok=True)
    if not sample_dirs:
        sys.stderr.write("postprocess: no sample directories found\n")

    records = []
    samples = []
    for sd in sample_dirs:
        name = args.sample if (args.sample_dir and args.sample) else \
            os.path.basename(os.path.normpath(sd))
        samples.append(name)
        recs = parsers.parse_sample_dir(sd, name)
        records.extend(recs)
        sys.stderr.write("postprocess: %s -> %d calls\n" % (name, len(recs)))

    write_calls(records, os.path.join(out_dir, "cohort_sv_calls.tsv"))
    write_common_deletion(records, samples, os.path.join(out_dir, "cohort_common_deletion.tsv"))
    write_matrix(records, samples, os.path.join(out_dir, "cohort_caller_matrix.tsv"))
    write_summary(records, samples, os.path.join(out_dir, "cohort_summary.txt"))

    sys.stderr.write("postprocess: wrote cohort_*.tsv to %s\n" % out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
