#!/usr/bin/env python3
"""Evaluate cohort SV calls against the MitoHPC truth scenarios.

Mirrors the scenario coverage of the MitoHPC `sv-calling` test suite (the source
of our mock BAMs): across the diverse constructs — common deletion at varying
VAF/depth, a non-repeat deletion, a D-loop deletion, a multi-deletion sample, a
tandem duplication, an origin-crossing deletion, and wild-type — it checks that
our 5-caller suite is

  * SENSITIVE   — the truth deletion is detected by >=1 caller, and
  * SPECIFIC    — no caller spuriously calls the ~4977 bp COMMON deletion on a
                  sample whose truth does not contain it (incl. the duplication,
                  origin-crossing and wild-type negatives).

Because our callers are heuristic (and quite different from MitoHPC's own
caller), only a small, reliable set of expectations are HARD gates; the rest are
reported as warnings in a scenario x caller matrix so coverage is visible without
making CI brittle.

Inputs:  cohort_sv_calls.tsv (from postprocess.py) + truth.tsv (+ optional real
truth). Output: a markdown scenario matrix, and a non-zero exit on any HARD-gate
failure.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "pipeline", "lib"))
import parsers  # noqa: E402

GEN_TOL = 250   # generic per-breakpoint match tolerance (bp), cross-caller

# Samples where the truth common-deletion is reliable enough across our callers
# to gate the build on (verified detected by eKLIPse/MitoMut/MitoSeek).
HARD_SENSITIVITY = {"sv_del4977_h30"}

# Real-data expectations (not in the mock truth.tsv); committed real BAMs.
REAL_TRUTH = {
    "spike_del4977_h20": [("del", 8469, 13447)],   # del4977 spiked @~20% -> sensitivity
    "NA12718":           [("none", None, None)],    # healthy 1000G -> specificity
    "NA12748":           [("none", None, None)],
    "NA12775":           [("none", None, None)],
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
                truth.setdefault(f[0], []).append((f[1], _num(f[2]), _num(f[3])))
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


def _match(c5, c3, e5, e3, tol=GEN_TOL):
    if None in (c5, c3, e5, e3):
        return False
    return abs(c5 - e5) <= tol and abs(c3 - e3) <= tol


def detectors(calls, sample, e5, e3, common):
    """Callers whose deletion call matches the expected event."""
    out = set()
    for c in calls.get(sample, []):
        if c["sv_type"] not in ("deletion", "duplication"):
            continue
        if common:
            if parsers.is_common_deletion(c["bp5"], c["bp3"]):
                out.add(c["caller"])
        elif _match(c["bp5"], c["bp3"], e5, e3):
            out.add(c["caller"])
    return out


def common_callers(calls, sample):
    return {c["caller"] for c in calls.get(sample, [])
            if parsers.is_common_deletion(c["bp5"], c["bp3"])}


def base_sample(s):
    """Map a CRAM round-trip sample back to its BAM truth (sv_x_cram -> sv_x)."""
    return s[:-5] if s.endswith("_cram") else s


def evaluate(truth, calls, run_samples=None):
    """Return (rows, hard_failures, warnings). rows = matrix lines.

    If run_samples is given, evaluate exactly those samples (so a 0-call
    wild-type sample, absent from cohort_sv_calls.tsv, is still scored and
    samples we never ran are not). Otherwise fall back to truth + cohort.
    """
    rows, hard, warn = [], [], []
    if run_samples:
        order = list(run_samples)
    else:
        order = list(truth.keys())
        for s in calls:
            if s not in order and base_sample(s) in truth:
                order.append(s)
    seen = []
    for s in order:
        if s and s not in seen:
            seen.append(s)

    for sample in seen:
        events = truth.get(base_sample(sample), truth.get(sample, []))
        carries_common = any(k in ("del", "delwrap")
                             and parsers.is_common_deletion(b5, b3)
                             for (k, b5, b3) in events)

        # --- per-event sensitivity ---
        if not events or all(e[0] == "none" for e in events):
            rows.append((sample, "wild-type (no SV)", "-",
                         "n/a (specificity sample)"))
        for (kind, e5, e3) in events:
            if kind in ("del", "delwrap"):
                is_common = parsers.is_common_deletion(e5, e3)
                det = detectors(calls, sample, e5, e3, common=is_common)
                label = "del %s-%s%s" % (e5, e3, " [COMMON]" if is_common else "")
                rows.append((sample, label, "yes" if det else "NO",
                             ", ".join(sorted(det)) if det else "(none)"))
                if not det:
                    msg = "sensitivity: %s del %s-%s not detected by any caller" % (
                        sample, e5, e3)
                    if base_sample(sample) in HARD_SENSITIVITY and is_common:
                        hard.append(msg)
                    else:
                        warn.append(msg)
            elif kind == "dup":
                det = detectors(calls, sample, e5, e3, common=False)
                rows.append((sample, "dup %s-%s" % (e5, e3),
                             "yes" if det else "no",
                             ", ".join(sorted(det)) if det else "(none)"))

        # --- specificity: no spurious COMMON deletion on non-carriers ---
        if not carries_common:
            fp = common_callers(calls, sample)
            if fp:
                hard.append("specificity: %s does not carry the common deletion "
                            "but it was called by %s" % (sample, ", ".join(sorted(fp))))
    return rows, hard, warn


def write_matrix(rows, hard, warn, path):
    with open(path, "w") as fh:
        fh.write("## Scenario coverage (vs MitoHPC truth)\n\n")
        fh.write("One row per truth event across the mock + real test cohort. "
                 "**detected** = at least one caller matched the truth deletion "
                 "(common deletion within +/-%d bp; others within +/-%d bp).\n\n"
                 % (parsers.COMMON_DEL_TOL, GEN_TOL))
        fh.write("| sample | truth event | detected | callers |\n")
        fh.write("|--------|-------------|:--------:|---------|\n")
        for (s, ev, det, who) in rows:
            fh.write("| %s | %s | %s | %s |\n" % (s, ev, det, who))
        fh.write("\n")
        fh.write("HARD-gate failures: %d  |  warnings: %d\n\n" % (len(hard), len(warn)))
        if hard:
            fh.write("**HARD failures:**\n\n")
            for h in hard:
                fh.write("- %s\n" % h)
            fh.write("\n")
        if warn:
            fh.write("**Warnings (heuristic-caller sensitivity misses, not gated):**\n\n")
            for w in warn:
                fh.write("- %s\n" % w)
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
    if not args.no_real:
        truth.update(REAL_TRUTH)
    calls = load_calls(args.calls)
    run_samples = [s for s in args.samples.replace(",", " ").split() if s] or None

    rows, hard, warn = evaluate(truth, calls, run_samples)
    write_matrix(rows, hard, warn, args.out_md)

    # Console report
    sys.stderr.write("\n=== scenario evaluation ===\n")
    for (s, ev, det, who) in rows:
        sys.stderr.write("  [%s] %-22s %-30s -> %s (%s)\n"
                         % ("OK " if det in ("yes", "-") else "   ", s, ev, det, who))
    for w in warn:
        sys.stderr.write("WARNING: %s\n" % w)
    for h in hard:
        sys.stderr.write("HARD-FAIL: %s\n" % h)
    sys.stderr.write("scenario gates: %d hard failure(s), %d warning(s)\n"
                     % (len(hard), len(warn)))
    return 1 if hard else 0


if __name__ == "__main__":
    raise SystemExit(main())
