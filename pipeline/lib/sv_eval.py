"""Shared SV-scenario evaluation: categorize truth, score detections, and compute
per-category / overall accuracy metrics for the caller comparison.

EVALUATION ONLY — these numbers never gate the build. We treat all six callers as
binary *deletion* detectors and score them against the MitoHPC labelled cohort.

Classes (per truth event for deletions; per sample for negatives):
  positive  — a real, detectable deletion event (kind=del, expect pass/detected,
              or the real del4977 spike). A matching deletion-like call => TP.
  negative  — NO deletion present, so a deletion-like call => false positive
              (kind none/dup/inv/invdup; the real 1000G healthy samples).
  excluded  — evaluation-only rows kept OUT of the confusion matrix:
                sub     kind=del expect=no_record (below detectable size),
                wrap    kind=delwrap (origin-crossing; reported as the complement),
                knownfp kind=dupdel (a hit is a documented gap, not credit).

A "deletion-like call" = a parsed record typed deletion/duplication with a span
(bp5 & bp3). duplication is included because some callers mis-type the common
deletion and because a span-call on a dup IS a false deletion for a del detector.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parsers  # noqa: E402

GEN_TOL = 250                    # generic per-breakpoint match tolerance (bp)

KIND2CAT = {"del": "del", "delwrap": "origin", "dup": "dup", "inv": "inv",
            "dupdel": "complex", "invdup": "complex", "none": "control"}
# Tab order; "all" is appended by the report. Deletions first = the default suite.
CATEGORY_ORDER = ["del", "control", "dup", "inv", "origin", "complex"]
CATEGORY_LABEL = {"del": "Deletions", "control": "Controls", "dup": "Duplications",
                  "inv": "Inversions", "origin": "Origin-crossing",
                  "complex": "Complex", "all": "All"}
# Categories that are forward-looking / ambiguous -> shown but flagged eval-only.
EVAL_ONLY_CATS = {"dup", "inv", "origin", "complex"}


def category_of(kind):
    return KIND2CAT.get(kind, "complex")


def classify(kind, expect):
    """(klass, reason) for a truth event/sample. klass in positive|negative|excluded."""
    if kind == "del":
        if expect == "no_record":
            return ("excluded", "sub")
        return ("positive", "")                       # pass/detected/'' (real spike)
    if kind == "delwrap":
        return ("excluded", "wrap")
    if kind == "dupdel":
        return ("excluded", "knownfp")
    if kind in ("dup", "inv", "invdup", "none"):
        return ("negative", "")
    return ("excluded", "other")


def _match(c5, c3, e5, e3, tol=GEN_TOL):
    if None in (c5, c3, e5, e3):
        return False
    return abs(c5 - e5) <= tol and abs(c3 - e3) <= tol


def deletion_like(call):
    return (call.get("sv_type") in ("deletion", "duplication")
            and call.get("bp5") is not None and call.get("bp3") is not None)


def build_trials(truth, samples, calls_by_sample, real_truth=None):
    """One trial per positive deletion event, per negative sample, and per excluded
    row. Each trial records, per caller, detection (positives/excluded) or a false
    positive (negatives). `truth` maps sample -> [(kind, bp5, bp3, expect?)].
    Tuples may be 3- or 4-long (expect optional)."""
    real_truth = real_truth or {}
    trials = []
    for sample in samples:
        base = sample[:-5] if sample.endswith("_cram") else sample
        events = truth.get(base, truth.get(sample, real_truth.get(base, [])))
        scalls = [c for c in calls_by_sample.get(sample, []) if deletion_like(c)]

        neg_kinds = {e[0] for e in events} if events else {"none"}
        # A sample with no del/delwrap event is a single NEGATIVE specificity trial.
        has_real_del = any(e[0] in ("del", "delwrap") for e in events)
        if not has_real_del:
            kind = events[0][0] if events else "none"
            expect = events[0][3] if (events and len(events[0]) > 3) else ""
            klass, reason = classify(kind, expect)
            fp = {}
            for c in scalls:
                common = parsers.is_common_deletion(c["bp5"], c["bp3"])
                fp[c["caller"]] = {"common": bool(common)}
            trials.append({"sample": sample, "category": category_of(kind),
                           "kind": kind, "klass": klass, "reason": reason,
                           "label": "wild-type (no SV)" if kind == "none"
                                    else "%s (no deletion)" % kind,
                           "is_common": False, "fp": fp, "detected": {},
                           "eval_only": category_of(kind) in EVAL_ONLY_CATS})
            continue

        # Deletion / delwrap events: one trial each, greedily matched.
        used = set()
        for ev in events:
            kind, e5, e3 = ev[0], ev[1], ev[2]
            expect = ev[3] if len(ev) > 3 else ""
            if kind not in ("del", "delwrap"):
                continue
            klass, reason = classify(kind, expect)
            is_common = parsers.is_common_deletion(e5, e3)
            det = {}
            for c in scalls:
                key = (c["caller"], c.get("bp5"), c.get("bp3"))
                if key in used:
                    continue
                hit = (parsers.is_common_deletion(c["bp5"], c["bp3"]) if is_common
                       else _match(c.get("bp5"), c.get("bp3"), e5, e3))
                if hit:
                    det[c["caller"]] = True
                    used.add(key)
            label = "%s %s–%s%s" % (kind, e5, e3, " · COMMON" if is_common else "")
            trials.append({"sample": sample, "category": category_of(kind),
                           "kind": kind, "klass": klass, "reason": reason,
                           "label": label, "is_common": bool(is_common),
                           "detected": det, "fp": {},
                           "eval_only": klass == "excluded"})
    return trials


def _metrics(tp, fn, fp, tn):
    import math
    rec = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    prec = tp / (tp + fp) if (tp + fp) else None
    f1 = (2 * prec * rec / (prec + rec)) if (prec and rec and (prec + rec)) else None
    bacc = ((rec + spec) / 2) if (rec is not None and spec is not None) else None
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = ((tp * tn - fp * fn) / math.sqrt(denom)) if denom else None
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "sensitivity": rec, "specificity": spec, "precision": prec,
            "f1": f1, "balanced_acc": bacc, "mcc": mcc}


def metrics_by_category(trials, callers):
    """{caller: {category|'all': metrics+confusion+fp_samples}}. Excluded trials
    are not scored. Negatives contribute TN/FP, positives TP/FN."""
    out = {c: {} for c in callers}
    cats = sorted({t["category"] for t in trials}) + ["all"]
    for c in callers:
        for cat in cats:
            tp = fn = fp = tn = 0
            fp_samples, common_fp = [], []
            for t in trials:
                if cat != "all" and t["category"] != cat:
                    continue
                if t["klass"] == "positive":
                    if t["detected"].get(c):
                        tp += 1
                    else:
                        fn += 1
                elif t["klass"] == "negative":
                    info = t["fp"].get(c)
                    if info:
                        fp += 1
                        fp_samples.append(t["sample"])
                        if info.get("common"):
                            common_fp.append(t["sample"])
                    else:
                        tn += 1
            m = _metrics(tp, fn, fp, tn)
            m["fp_samples"] = sorted(set(fp_samples))
            m["common_fp_samples"] = sorted(set(common_fp))
            m["n_pos"] = tp + fn
            m["n_neg"] = tn + fp
            out[c][cat] = m
    return out
