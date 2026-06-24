#!/usr/bin/env python3
"""Unit tests for the scenario evaluator (test/check_scenarios.py).

Runs without Docker: builds synthetic truth + cohort-call tables and asserts the
HARD specificity/sensitivity gates and the per-scenario warnings behave.

Run:  python3 test/test_check_scenarios.py   (or: pytest test/test_check_scenarios.py)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_scenarios as cs  # noqa: E402

TRUTH = {
    "sv_del4977_h30": [("del", 8469, 13447)],          # common deletion
    "sv_del6000_h50": [("del", 5999, 10999)],          # non-repeat deletion
    "sv_dup":         [("dup", 6000, 7000)],           # duplication (no common del)
    "sv_origin":      [("delwrap", 16400, 200)],       # origin-crossing
    "sv_wt":          [("none", None, None)],          # wild-type
}


def _calls(rows):
    out = {}
    for sample, caller, sv_type, bp5, bp3 in rows:
        out.setdefault(sample, []).append(
            {"caller": caller, "sv_type": sv_type, "bp5": bp5, "bp3": bp3})
    return out


def test_good_cohort_passes():
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_del4977_h30", "mitomut", "deletion", 8482, 13446),
        ("sv_del6000_h50", "eklipse", "deletion", 5999, 10999),
        ("sv_dup", "mitoseek", "deletion", 6100, 6960),     # a non-common call: fine
    ])
    samples = ["sv_del4977_h30", "sv_del6000_h50", "sv_dup", "sv_origin", "sv_wt"]
    rows, hard, warn = cs.evaluate(TRUTH, calls, samples)
    assert hard == [], hard
    # origin not detected -> a (non-gated) warning
    assert any("sv_origin" in w for w in warn)


def test_specificity_hard_fail_on_wildtype_common_call():
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_wt", "eklipse", "deletion", 8470, 13447),       # spurious COMMON on WT
    ])
    rows, hard, warn = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert any("specificity" in h and "sv_wt" in h for h in hard), hard


def test_specificity_hard_fail_on_duplication_common_call():
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_dup", "mitomut", "deletion", 8469, 13447),      # spurious COMMON on dup
    ])
    rows, hard, warn = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert any("specificity" in h and "sv_dup" in h for h in hard), hard


def test_sensitivity_hard_fail_when_canonical_missing():
    calls = _calls([
        ("sv_del6000_h50", "eklipse", "deletion", 5999, 10999),
    ])  # nothing for sv_del4977_h30
    rows, hard, warn = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert any("sensitivity" in h and "sv_del4977_h30" in h for h in hard), hard


def test_non_common_deletion_does_not_trip_specificity():
    # del6000 carries a (non-common) deletion; detecting it must NOT be a
    # specificity failure, and not detecting the common deletion is fine.
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_del6000_h50", "eklipse", "deletion", 5999, 10999),
    ])
    rows, hard, warn = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert hard == [], hard


def test_cram_alias_inherits_truth():
    truth = {"sv_del4977_h30": [("del", 8469, 13447)]}
    calls = _calls([("sv_del4977_h30_cram", "mitomut", "deletion", 8482, 13446)])
    rows, hard, warn = cs.evaluate(truth, calls, ["sv_del4977_h30_cram"])
    # the cram sample is scored against the base sample's common-deletion truth
    assert any(r[0] == "sv_del4977_h30_cram" and r[2] == "yes" for r in rows), rows


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
