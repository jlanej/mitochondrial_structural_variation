#!/usr/bin/env python3
"""Unit tests for the scenario comparison tool (test/check_scenarios.py).

The tool is EVALUATION-ONLY: it classifies, per scenario, which callers detect
each truth deletion and where a caller reports the common deletion on a
non-carrier — but it never gates (always exits 0). These tests assert the
classification (sensitivity misses / specificity observations) and that it does
not fail the build.

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


def test_detection_recorded_and_nothing_gated():
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_del4977_h30", "mitomut", "deletion", 8482, 13446),
        ("sv_del6000_h50", "eklipse", "deletion", 5999, 10999),
    ])
    rows, sens_miss, spec_obs = cs.evaluate(TRUTH, calls, list(TRUTH))
    # del4977_h30 detected by two callers -> a "yes" row listing both
    assert any(r[0] == "sv_del4977_h30" and r[2] == "yes"
               and "eklipse" in r[3] and "mitomut" in r[3] for r in rows), rows
    # origin not detected -> recorded as a sensitivity observation (not a gate)
    assert any("sv_origin" in m for m in sens_miss)
    # no false-positive common-deletion calls
    assert spec_obs == [], spec_obs


def test_wildtype_common_call_is_observation_only():
    calls = _calls([
        ("sv_del4977_h30", "eklipse", "deletion", 8469, 13447),
        ("sv_wt", "eklipse", "deletion", 8470, 13447),       # FP common del on WT
    ])
    rows, sens_miss, spec_obs = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert any("sv_wt" in m and "eklipse" in m for m in spec_obs), spec_obs


def test_duplication_common_call_is_observation_only():
    calls = _calls([
        ("sv_dup", "mitomut", "deletion", 8469, 13447),      # FP common del on dup
    ])
    rows, sens_miss, spec_obs = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert any("sv_dup" in m for m in spec_obs), spec_obs


def test_non_common_deletion_is_not_a_specificity_observation():
    # del6000 carries a (non-common) deletion; detecting it must NOT be flagged
    # as a common-deletion false positive.
    calls = _calls([("sv_del6000_h50", "eklipse", "deletion", 5999, 10999)])
    rows, sens_miss, spec_obs = cs.evaluate(TRUTH, calls, list(TRUTH))
    assert spec_obs == [], spec_obs


def test_cram_alias_inherits_truth():
    truth = {"sv_del4977_h30": [("del", 8469, 13447)]}
    calls = _calls([("sv_del4977_h30_cram", "mitomut", "deletion", 8482, 13446)])
    rows, sens_miss, spec_obs = cs.evaluate(truth, calls, ["sv_del4977_h30_cram"])
    assert any(r[0] == "sv_del4977_h30_cram" and r[2] == "yes" for r in rows), rows


def test_main_always_exits_zero(tmp_path=None):
    import tempfile
    d = tempfile.mkdtemp()
    calls = os.path.join(d, "calls.tsv")
    truth = os.path.join(d, "truth.tsv")
    out = os.path.join(d, "m.md")
    with open(calls, "w") as fh:
        fh.write("sample\tcaller\tsv_type\tbp5\tbp3\tsvlen\tsupport\thet\tcommon_deletion\textra\n")
        # a wild-type false positive — would be a "gate" failure if we gated, but
        # this tool only reports, so it must still exit 0.
        fh.write("sv_wt\teklipse\tdeletion\t8470\t13447\t4977\t5\t0.01\t1\tx\n")
    with open(truth, "w") as fh:
        fh.write("#sample\tkind\tbp5\tbp3\tsvlen\thet\tdepth\n")
        fh.write("sv_wt\tnone\t.\t.\t.\t0\t300\n")
    rc = cs.main(["--calls", calls, "--truth", truth, "--out-md", out,
                  "--samples", "sv_wt", "--no-real"])
    assert rc == 0, rc
    assert os.path.isfile(out)
    import shutil
    shutil.rmtree(d, ignore_errors=True)


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
