#!/usr/bin/env python3
"""Unit tests for the shared SV-scenario evaluation kernel (pipeline/lib/sv_eval.py).

Covers categorization, trial classification (positive / negative / excluded),
detection + false-positive scoring, the multidel greedy match, the CRAM alias,
and the per-category / overall accuracy metrics.

Run:  python3 test/test_sv_eval.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lib"))
import sv_eval as E  # noqa: E402

CALLERS = ["mitohpc", "eklipse", "mitomut", "mitoseek"]


def call(c, t="deletion", b5=None, b3=None):
    return {"caller": c, "sv_type": t, "bp5": b5, "bp3": b3}


def trials_for(truth, calls_by_sample, samples=None, real=None):
    samples = samples or list(truth) + list(real or {})
    return E.build_trials(truth, samples, calls_by_sample, real or {})


def test_category_and_classify():
    assert E.category_of("del") == "del" and E.category_of("delwrap") == "origin"
    assert E.category_of("dup") == "dup" and E.category_of("invdup") == "complex"
    assert E.classify("del", "pass") == ("positive", "")
    assert E.classify("del", "no_record") == ("excluded", "sub")
    assert E.classify("delwrap", "wrap") == ("excluded", "wrap")
    assert E.classify("dupdel", "known_fp") == ("excluded", "knownfp")
    assert E.classify("dup", "no_pass")[0] == "negative"
    assert E.classify("none", "no_pass")[0] == "negative"


def test_positive_detection_and_miss():
    truth = {"sv_del4977_h30": [("del", 8469, 13447, "pass")]}
    cbs = {"sv_del4977_h30": [call("eklipse", b5=8469, b3=13447),
                              call("mitomut", b5=8482, b3=13446)]}
    t = trials_for(truth, cbs)[0]
    assert t["category"] == "del" and t["klass"] == "positive" and t["is_common"]
    assert t["detected"].get("eklipse") and t["detected"].get("mitomut")
    assert "mitoseek" not in t["detected"]


def test_false_positive_on_negatives():
    truth = {"sv_wt": [("none", None, None, "no_pass")],
             "sv_dup": [("dup", 6000, 7000, "no_pass")]}
    cbs = {"sv_wt": [call("eklipse", b5=8470, b3=13447)],          # common-del FP
           "sv_dup": [call("mitomut", b5=6000, b3=7000)]}          # del call on a dup
    ts = {t["sample"]: t for t in trials_for(truth, cbs)}
    assert ts["sv_wt"]["klass"] == "negative"
    assert ts["sv_wt"]["fp"]["eklipse"]["common"] is True
    assert ts["sv_dup"]["fp"]["mitomut"] and not ts["sv_dup"]["fp"]["mitomut"]["common"]


def test_excluded_rows_not_scored():
    truth = {"sv_origin": [("delwrap", 16400, 200, "wrap")],
             "sv_del_45": [("del", 9000, 9046, "no_record")],
             "sv_dupdel": [("dupdel", 5000, 8000, "known_fp")]}
    ts = {t["sample"]: t for t in trials_for(truth, {})}
    assert ts["sv_origin"]["klass"] == "excluded" and ts["sv_origin"]["reason"] == "wrap"
    assert ts["sv_del_45"]["reason"] == "sub" and ts["sv_del_45"]["eval_only"]
    assert ts["sv_dupdel"]["reason"] == "knownfp" and ts["sv_dupdel"]["category"] == "complex"
    m = E.metrics_by_category(list(ts.values()), CALLERS)
    # nothing scored -> all confusion counts zero
    assert m["eklipse"]["all"]["tp"] == 0 and m["eklipse"]["all"]["n_neg"] == 0


def test_multidel_two_events_greedy():
    truth = {"sv_multidel": [("del", 8469, 13447, "pass"), ("del", 5999, 10999, "pass")]}
    # one caller finds both (distinct records), another only del4977
    cbs = {"sv_multidel": [call("mitohpc", b5=8482, b3=13446),
                           call("mitohpc", b5=6000, b3=10999),
                           call("eklipse", b5=8470, b3=13447)]}
    ts = trials_for(truth, cbs)
    assert len(ts) == 2
    det = [t["detected"] for t in ts]
    assert det[0].get("mitohpc") and det[1].get("mitohpc")     # both events
    # eklipse matched at most ONE event (the common one), not both
    assert sum(1 for d in det if d.get("eklipse")) == 1


def test_cram_alias_inherits_truth():
    truth = {"sv_del4977_h30": [("del", 8469, 13447, "pass")]}
    cbs = {"sv_del4977_h30_cram": [call("mitomut", b5=8482, b3=13446)]}
    t = trials_for(truth, cbs, samples=["sv_del4977_h30_cram"])[0]
    assert t["sample"] == "sv_del4977_h30_cram" and t["detected"].get("mitomut")


def test_non_deletion_call_ignored():
    # a mitoseek 'breakpoint' row (no span) is neither a detection nor an FP
    truth = {"sv_wt": [("none", None, None, "no_pass")]}
    cbs = {"sv_wt": [call("mitoseek", t="breakpoint", b5=100, b3=None)]}
    t = trials_for(truth, cbs)[0]
    assert "mitoseek" not in t["fp"]


def test_metrics_confusion_and_mcc():
    truth = {
        "sv_del4977_h30": [("del", 8469, 13447, "pass")],
        "sv_del6000_h50": [("del", 5999, 10999, "pass")],
        "sv_wt": [("none", None, None, "no_pass")],
        "sv_dup": [("dup", 6000, 7000, "no_pass")],
    }
    cbs = {
        "sv_del4977_h30": [call("mitohpc", b5=8482, b3=13446)],   # TP
        "sv_del6000_h50": [call("mitohpc", b5=5999, b3=10999)],   # TP
        "sv_wt": [call("mitohpc", b5=1000, b3=2000)],             # FP
    }
    m = E.metrics_by_category(trials_for(truth, cbs), CALLERS)["mitohpc"]
    a = m["all"]
    assert (a["tp"], a["fn"], a["fp"], a["tn"]) == (2, 0, 1, 1)
    assert abs(a["sensitivity"] - 1.0) < 1e-9
    assert abs(a["specificity"] - 0.5) < 1e-9
    assert a["mcc"] is not None and math.isfinite(a["mcc"])
    # per-category: deletions has only positives (no specificity), controls only neg
    assert m["del"]["n_pos"] == 2 and m["del"]["n_neg"] == 0
    assert m["control"]["n_neg"] == 1 and m["control"]["specificity"] == 0.0
    assert m["dup"]["specificity"] == 1.0          # no FP on the dup


def _main():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn(); print("PASS", fn.__name__)
        except AssertionError as e:
            failed += 1; print("FAIL", fn.__name__, "->", e)
        except Exception as e:  # noqa: BLE001
            failed += 1; print("ERROR", fn.__name__, "->", repr(e))
    print("\n%d/%d passed" % (len(fns) - failed, len(fns)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())
