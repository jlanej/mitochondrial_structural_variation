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

GEN_TOL = 250   # generic per-breakpoint match tolerance (bp), cross-caller

# Real-data expectations (not in the mock truth.tsv); committed real BAMs.
REAL_TRUTH = {
    "spike_del4977_h20": [("del", 8469, 13447)],   # del4977 spiked @~20%
    "NA12718":           [("none", None, None)],    # healthy 1000G
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
    """Callers whose deletion/duplication call matches the expected event."""
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
    """Return (rows, sensitivity_misses, specificity_observations).

    All three are INFORMATIONAL — nothing here gates the build. If run_samples is
    given, score exactly those samples (so a 0-call wild-type sample, absent from
    cohort_sv_calls.tsv, is still scored and samples we never ran are not).
    """
    rows, sens_miss, spec_obs = [], [], []
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

        if not events or all(e[0] == "none" for e in events):
            rows.append((sample, "wild-type (no SV)", "-",
                         "n/a (specificity sample)"))
        for (kind, e5, e3) in events:
            if kind in ("del", "delwrap"):
                is_common = parsers.is_common_deletion(e5, e3)
                det = detectors(calls, sample, e5, e3, common=is_common)
                label = "del %s-%s%s" % (e5, e3, " [COMMON]" if is_common else "")
                rows.append((sample, label, "yes" if det else "no",
                             ", ".join(sorted(det)) if det else "(none)"))
                if not det:
                    sens_miss.append("%s del %s-%s detected by no caller" % (sample, e5, e3))
            elif kind == "dup":
                det = detectors(calls, sample, e5, e3, common=False)
                rows.append((sample, "dup %s-%s" % (e5, e3),
                             "yes" if det else "no",
                             ", ".join(sorted(det)) if det else "(none)"))

        # Observation only: a caller reporting the common deletion on a sample
        # that does not carry it (false-positive behaviour, recorded not gated).
        if not carries_common:
            fp = common_callers(calls, sample)
            if fp:
                spec_obs.append("%s does not carry the common deletion, "
                                "but it was reported by %s" % (sample, ", ".join(sorted(fp))))
    return rows, sens_miss, spec_obs


def write_matrix(rows, sens_miss, spec_obs, path):
    with open(path, "w") as fh:
        fh.write("## Caller comparison across MitoHPC scenarios\n\n")
        fh.write("Evaluation only — a record of how each third-party caller behaves "
                 "on the diverse test constructs (we do not control their source, "
                 "so nothing here gates the build). **detected** = at least one "
                 "caller matched the truth deletion (common deletion within "
                 "+/-%d bp; others within +/-%d bp).\n\n"
                 % (parsers.COMMON_DEL_TOL, GEN_TOL))
        fh.write("| sample | truth event | detected | callers |\n")
        fh.write("|--------|-------------|:--------:|---------|\n")
        for (s, ev, det, who) in rows:
            fh.write("| %s | %s | %s | %s |\n" % (s, ev, det, who))
        fh.write("\n")
        if sens_miss:
            fh.write("**Sensitivity — truth events no caller detected (observation):**\n\n")
            for m in sens_miss:
                fh.write("- %s\n" % m)
            fh.write("\n")
        if spec_obs:
            fh.write("**Specificity — common-deletion calls on non-carriers (observation):**\n\n")
            for m in spec_obs:
                fh.write("- %s\n" % m)
            fh.write("\n")
        if not sens_miss and not spec_obs:
            fh.write("_Every truth deletion was detected by >=1 caller and no "
                     "common-deletion false positives were observed._\n\n")


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

    rows, sens_miss, spec_obs = evaluate(truth, calls, run_samples)
    write_matrix(rows, sens_miss, spec_obs, args.out_md)

    sys.stderr.write("\n=== caller comparison (evaluation only, not gated) ===\n")
    for (s, ev, det, who) in rows:
        sys.stderr.write("  %-22s %-30s detected=%-3s  %s\n" % (s, ev, det, who))
    for m in sens_miss:
        sys.stderr.write("  note (sensitivity): %s\n" % m)
    for m in spec_obs:
        sys.stderr.write("  note (specificity): %s\n" % m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
