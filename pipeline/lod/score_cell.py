#!/usr/bin/env python3
"""Score one LOD cell: for each caller, did it detect the truth deletion?

Reads a run_sample.sh output dir (per-caller calls) + the cell truth.tsv, and
emits one TSV row per caller with MitoHPC's detection rule:

  detected = a deletion call within BP_TOL=30 summed breakpoint error
             (|called_bp5 - truth_bp5| + |called_bp3 - truth_bp3|).

The 30-bp summed tolerance absorbs the repeat-mediated breakpoint shift and the
bp3 / bp3-1 (retained vs last-deleted) convention difference between callers, so
all six callers are judged identically and fairly. Runtime comes from status.tsv;
'passed' (MitoHPC FILTER==PASS) is recorded for MitoHPC only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import parsers  # noqa: E402

BP_TOL = 30
COLUMNS = ["arm", "caller", "variant", "vaf", "depth", "rep", "seed",
           "truth_bp5", "truth_bp3", "truth_svlen", "detected", "n_calls",
           "called_bp5", "called_bp3", "bp5_err", "bp3_err", "bp_err_abs",
           "svlen_err", "af", "support", "passed", "status", "runtime_s"]


def load_truth(path):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t",
                                  fieldnames=None):
            return {k.lstrip("#"): v for k, v in row.items()}
    return {}


def load_status(sample_dir):
    out = {}
    p = os.path.join(sample_dir, "status.tsv")
    if os.path.isfile(p):
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[r["caller"]] = (r.get("status", ""), r.get("seconds", ""))
    return out


def best_match(records, t5, t3):
    """Closest deletion call to truth by summed breakpoint error."""
    best = None
    for r in records:
        if r["sv_type"] not in ("deletion",):
            continue
        if r["bp5"] is None or r["bp3"] is None:
            continue
        d = abs(r["bp5"] - t5) + abs(r["bp3"] - t3)
        if best is None or d < best[0]:
            best = (d, r)
    return best  # (summed_err, record) or None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-dir", required=True, help="run_sample.sh output dir")
    ap.add_argument("--truth", required=True, help="cell truth.tsv from gen_cell")
    ap.add_argument("--arm", required=True, help="pipeline | circular")
    ap.add_argument("--out", required=True, help="TSV to append rows to")
    ap.add_argument("--sample", default=None)
    args = ap.parse_args(argv)

    t = load_truth(args.truth)
    t5, t3 = int(t["bp5"]), int(t["bp3"])
    tsvlen = int(t["svlen"])
    recs = parsers.parse_sample_dir(args.sample_dir, args.sample or t["variant"])
    by_caller = {}
    for r in recs:
        by_caller.setdefault(r["caller"], []).append(r)
    status = load_status(args.sample_dir)

    # Only score callers that ACTUALLY RAN this cell. status.tsv has one row per
    # requested caller (run_sample.sh), so it is the authoritative "ran" set. This
    # matters under the per-caller LOD split (run_sample is invoked with a single
    # --callers): without it, score_cell would emit a detected=0 row for every
    # caller that did NOT run, and the merged sweep would carry 5 bogus zero rows
    # per cell per arm — collapsing detection rates and corrupting the LOD fits.
    ran = [c for c in parsers.CALLERS if c in status or c in by_caller]

    new = not os.path.exists(args.out) or os.path.getsize(args.out) == 0
    with open(args.out, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if new:
            w.writerow(COLUMNS)
        for caller in ran:
            crecs = by_caller.get(caller, [])
            m = best_match(crecs, t5, t3)
            st, secs = status.get(caller, ("", ""))
            if m is not None:
                d, r = m
                detected = 1 if d <= BP_TOL else 0
                passed = ("1" if "filter=PASS" in (r.get("extra") or "") else "0") \
                    if caller == "mitohpc" else "NA"
                w.writerow([args.arm, caller, t["variant"], t["vaf"], t["depth"],
                            t["rep"], t["seed"], t5, t3, tsvlen, detected, len(crecs),
                            r["bp5"], r["bp3"], r["bp5"] - t5, r["bp3"] - t3, d,
                            (r["svlen"] - tsvlen) if r["svlen"] is not None else "",
                            "" if r["het"] is None else round(r["het"], 4),
                            "" if r["support"] is None else r["support"],
                            passed, st, secs])
            else:
                w.writerow([args.arm, caller, t["variant"], t["vaf"], t["depth"],
                            t["rep"], t["seed"], t5, t3, tsvlen, 0, len(crecs),
                            "", "", "", "", "", "", "", "",
                            ("0" if caller == "mitohpc" else "NA"), st, secs])
    sys.stderr.write("[score_cell] arm=%s -> %d caller rows (%s)\n"
                     % (args.arm, len(ran), ",".join(ran) or "none"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
