#!/usr/bin/env python3
"""Compare what each caller does across the MitoHPC truth scenarios.

This is an EVALUATION / COMPARISON tool, not a pass/fail gate. We do not control
the third-party callers' source, so we do not assert that they must (or must
not) call any particular event — we simply record, for the diverse constructs in
the MitoHPC `sv-calling` test cohort (common deletion at varying VAF/depth, a
non-repeat deletion, a D-loop deletion, a multi-deletion, a tandem duplication,
an origin-crossing deletion, wild-type, plus real 1000G + a del4977 spike-in),
which callers detect each truth deletion and where a caller reports the ~4977 bp
common deletion on a sample that does not carry it.

The result is a scenario x caller matrix (markdown) plus two informational note
lists (sensitivity misses, specificity observations). It always exits 0 — the
build is gated elsewhere (callers run + produce output, post-processing works,
degenerate inputs fail cleanly), on things the pipeline controls.

Inputs: cohort_sv_calls.tsv (from postprocess.py) + truth.tsv (+ committed real
expectations). Output: the markdown matrix.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "pipeline", "lib"))
import parsers  # noqa: E402
import sv_eval  # noqa: E402

GEN_TOL = 250   # generic per-breakpoint match tolerance (bp), cross-caller

# Real-data expectations (not in the mock truth.tsv); committed real BAMs.
REAL_TRUTH = {
    "spike_del4977_h20": [("del", 8469, 13447, "pass")],   # del4977 spiked @~20%
    "NA12718":           [("none", None, None, "no_pass")],  # healthy 1000G
    "NA12748":           [("none", None, None, "no_pass")],
    "NA12775":           [("none", None, None, "no_pass")],
}


def _num(x):
    try:
        return int(round(float(x)))
    except (ValueError, TypeError):
        return None


def load_truth(path):
    truth = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                expect = f[7] if len(f) > 7 else ""
                truth.setdefault(f[0], []).append((f[1], _num(f[2]), _num(f[3]), expect))
    return truth


def load_calls(path):
    calls = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                calls.setdefault(row["sample"], []).append({
                    "caller": row["caller"], "sv_type": row["sv_type"],
                    "bp5": _num(row["bp5"]), "bp3": _num(row["bp3"]),
                })
    return calls


def _pct(v):
    return "—" if v is None else "%d%%" % round(100 * v)


def write_matrix(trials, metrics, callers, path):
    """Per-category scenario x caller matrix + per-caller accuracy table.
    EVALUATION ONLY — nothing here gates the build."""
    cats = [c for c in sv_eval.CATEGORY_ORDER
            if any(t["category"] == c for t in trials)]
    with open(path, "w") as fh:
        fh.write("## Caller comparison across MitoHPC scenarios\n\n")
        fh.write("Evaluation only — how each third-party caller behaves on the diverse "
                 "MitoHPC test constructs (we do not control their source, so nothing "
                 "here gates the build). A **deletion-like call** = a record typed "
                 "deletion/duplication with a span matching the truth event "
                 "(common deletion within +/-%d bp; others within +/-%d bp). "
                 "Forward-looking / ambiguous rows (origin-crossing, dup-del, sub-size) "
                 "are shown but not scored.\n\n"
                 % (parsers.COMMON_DEL_TOL, GEN_TOL))

        # Overall accuracy table.
        fh.write("### Accuracy (all scored scenarios)\n\n")
        fh.write("| caller | sensitivity | specificity | precision | F1 | bal.acc | MCC | FP |\n")
        fh.write("|--------|:-----------:|:-----------:|:---------:|:--:|:-------:|:---:|:--:|\n")
        for c in sorted(callers, key=lambda c: -(metrics[c]["all"].get("balanced_acc") or -1)):
            m = metrics[c]["all"]
            f1 = "—" if m["f1"] is None else "%.2f" % m["f1"]
            mcc = "—" if m["mcc"] is None else "%.2f" % m["mcc"]
            fh.write("| %s | %s (%d/%d) | %s (%d/%d) | %s | %s | %s | %s | %d |\n" % (
                c, _pct(m["sensitivity"]), m["tp"], m["n_pos"],
                _pct(m["specificity"]), m["tn"], m["n_neg"],
                _pct(m["precision"]), f1, _pct(m["balanced_acc"]), mcc, m["fp"]))
        fh.write("\n")

        # Per-category matrices.
        for cat in cats:
            crows = [t for t in trials if t["category"] == cat]
            if not crows:
                continue
            label = sv_eval.CATEGORY_LABEL.get(cat, cat)
            note = " *(evaluation-only / forward-looking)*" if cat in sv_eval.EVAL_ONLY_CATS else ""
            fh.write("### %s%s\n\n" % (label, note))
            fh.write("| sample | truth event | %s |\n" % " | ".join(callers))
            fh.write("|--------|-------------|%s\n" % ("--|" * len(callers)))
            for t in crows:
                cells = []
                for c in callers:
                    if t["klass"] == "positive":
                        cells.append("detected" if t["detected"].get(c) else "·")
                    elif t["klass"] == "negative":
                        f = t["fp"].get(c)
                        cells.append(("FP*" if f.get("common") else "FP") if f else "·")
                    else:
                        called = t["detected"].get(c) or t["fp"].get(c)
                        cells.append(t["reason"] if called else "·")
                fh.write("| %s | %s | %s |\n" % (t["sample"], t["label"], " | ".join(cells)))
            fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="cohort_sv_calls.tsv")
    ap.add_argument("--truth", required=True, help="mock truth.tsv")
    ap.add_argument("--out-md", required=True, help="scenario matrix markdown out")
    ap.add_argument("--samples", default="",
                    help="space/comma list of samples actually run (restricts scoring)")
    ap.add_argument("--no-real", action="store_true",
                    help="do not add the committed real-BAM expectations")
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    real = {} if args.no_real else REAL_TRUTH
    calls = load_calls(args.calls)
    run_samples = [s for s in args.samples.replace(",", " ").split() if s]
    if not run_samples:
        run_samples = list(truth.keys()) + [s for s in real if s not in truth]
        for s in calls:
            base = s[:-5] if s.endswith("_cram") else s
            if s not in run_samples and base in truth:
                run_samples.append(s)
    callers = list(parsers.CALLERS)
    cbs = {s: list(cl) for s, cl in calls.items()}

    trials = sv_eval.build_trials(truth, run_samples, cbs, real)
    metrics = sv_eval.metrics_by_category(trials, callers)
    write_matrix(trials, metrics, callers, args.out_md)

    sys.stderr.write("\n=== caller comparison (evaluation only, not gated) ===\n")
    for c in callers:
        m = metrics[c]["all"]
        sys.stderr.write("  %-13s sens=%s (%d/%d)  spec=%s (%d/%d)  FP=%d  MCC=%s\n" % (
            c, _pct(m["sensitivity"]), m["tp"], m["n_pos"],
            _pct(m["specificity"]), m["tn"], m["n_neg"], m["fp"],
            "—" if m["mcc"] is None else "%.2f" % m["mcc"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
