#!/usr/bin/env python3
"""Aggregate the LOD sweep shards into per-cell rates and per-(caller,depth) LODs.

Input: the concatenated sweep TSV (score_cell schema; one row per
arm x caller x variant x vaf x depth x replicate). Output:

  lod_cells.tsv  one row per (arm,caller,variant,depth,vaf): detection k/n/rate
                 + Wilson 95% CI, MitoHPC PASS-rate, mean runtime.
  lod_fits.tsv   one row per (arm,caller,variant,depth): empirical transition +
                 reliable-detection heteroplasmy (headline), Firth-logistic
                 LOD50/LOD95 + bootstrap CI (supporting), near-separable flag.
  lod_runtime.tsv one row per (arm,caller) and per (all,caller): per-cell runtime
                 distribution (n, mean, median, p25/p75, min/max, total seconds).

Pure stdlib (uses lod_stats.py).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lod_stats as S  # noqa: E402

CELL_COLS = ["arm", "caller", "variant", "depth", "vaf", "det_k", "det_n",
             "det_rate", "det_lo", "det_hi", "pass_k", "pass_rate", "mean_runtime_s"]
FIT_COLS = ["arm", "caller", "variant", "depth", "n_levels", "n_reps",
            "emp_transition", "emp_reliable", "near_separable",
            "lod50", "lod95", "lod95_lo", "lod95_hi"]
RT_COLS = ["arm", "caller", "n", "mean_s", "median_s", "p25_s", "p75_s",
           "min_s", "max_s", "total_s"]


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def load_sweep(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def _fmt(v):
    return "" if v is None or (isinstance(v, float) and v != v) else (
        "%.4f" % v if isinstance(v, float) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", required=True, help="concatenated sweep TSV")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    rows = load_sweep(args.sweep)

    # group by (arm,caller,variant,depth,vaf) -> detected list, passed list, runtimes
    cells = defaultdict(lambda: {"det": [], "pass": [], "rt": []})
    for r in rows:
        key = (r["arm"], r["caller"], r["variant"], int(r["depth"]), float(r["vaf"]))
        c = cells[key]
        c["det"].append(1 if r.get("detected") == "1" else 0)
        if r.get("passed") not in ("NA", "", None):
            c["pass"].append(1 if r.get("passed") == "1" else 0)
        rt = _f(r.get("runtime_s"))
        if rt is not None and r.get("status") == "ok":
            c["rt"].append(rt)

    with open(os.path.join(args.outdir, "lod_cells.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(CELL_COLS)
        for (arm, caller, variant, depth, vaf) in sorted(cells):
            c = cells[(arm, caller, variant, depth, vaf)]
            k, n = sum(c["det"]), len(c["det"])
            p, lo, hi = S.wilson(k, n)
            pk, pn = sum(c["pass"]), len(c["pass"])
            prate = (pk / pn) if pn else None
            mrt = (sum(c["rt"]) / len(c["rt"])) if c["rt"] else None
            w.writerow([arm, caller, variant, depth, _fmt(vaf), k, n, _fmt(p),
                        _fmt(lo), _fmt(hi), pk, _fmt(prate), _fmt(mrt)])

    # group by (arm,caller,variant,depth) -> (vaf, detected) units for fitting
    groups = defaultdict(list)
    for r in rows:
        groups[(r["arm"], r["caller"], r["variant"], int(r["depth"]))].append(
            (float(r["vaf"]), 1 if r.get("detected") == "1" else 0))

    with open(os.path.join(args.outdir, "lod_fits.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(FIT_COLS)
        for (arm, caller, variant, depth) in sorted(groups):
            units = groups[(arm, caller, variant, depth)]
            xs = [u[0] for u in units]
            ys = [u[1] for u in units]
            # per-level rates for empirical readout
            lvl = defaultdict(list)
            for v, y in units:
                lvl[v].append(y)
            levels = [(v, sum(ys_) / len(ys_)) for v, ys_ in lvl.items()]
            counts = [(sum(ys_), len(ys_)) for _, ys_ in lvl.items()]
            emp = S.empirical_lod(levels)
            near_sep = S.empirical_separable(counts)
            beta = S.fit_logistic_firth([x for x in xs if True], ys)
            lod50 = S.lod_at(beta, 0.5)
            lod95 = S.lod_at(beta, 0.95)
            lo95, hi95 = S.cluster_bootstrap_lod(units, 0.95) if not near_sep else (S.NAN, S.NAN)
            w.writerow([arm, caller, variant, depth, len(lvl), len(units),
                        _fmt(emp["transition_hi"]), _fmt(emp["reliable_lo"]),
                        1 if near_sep else 0,
                        _fmt(lod50), _fmt(lod95), _fmt(lo95), _fmt(hi95)])

    # ---- runtime summary per (arm,caller) + per (all,caller) ----------------
    rt_by = defaultdict(list)          # (arm,caller) -> [runtime_s]
    rt_all = defaultdict(list)         # caller -> [runtime_s] across arms
    for r in rows:
        rt = _f(r.get("runtime_s"))
        if rt is None or r.get("status") != "ok":
            continue
        rt_by[(r["arm"], r["caller"])].append(rt)
        rt_all[r["caller"]].append(rt)

    def _rt_row(arm, caller, vals):
        b = S.boxstats(vals)
        if not b:
            return None
        return [arm, caller, b["n"], _fmt(b["mean"]), _fmt(b["med"]),
                _fmt(b["q1"]), _fmt(b["q3"]), _fmt(b["min"]), _fmt(b["max"]),
                _fmt(sum(vals))]

    with open(os.path.join(args.outdir, "lod_runtime.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(RT_COLS)
        for caller in sorted(rt_all):
            row = _rt_row("all", caller, rt_all[caller])
            if row:
                w.writerow(row)
        for (arm, caller) in sorted(rt_by):
            row = _rt_row(arm, caller, rt_by[(arm, caller)])
            if row:
                w.writerow(row)

    sys.stderr.write("[analyze_lod] %d cells, %d (caller,depth) fits, "
                     "%d runtime groups -> %s\n"
                     % (len(cells), len(groups), len(rt_all) + len(rt_by), args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
