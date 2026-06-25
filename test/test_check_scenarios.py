#!/usr/bin/env python3
"""Unit tests for the scenario comparison tool (test/check_scenarios.py).

The tool is EVALUATION-ONLY: it categorizes scenarios, scores each caller as a
deletion detector (per-category + overall accuracy), and renders a markdown
matrix — but it never gates (always exits 0). The detailed scoring logic is
covered in test_sv_eval.py; here we assert the tool's integration + that it does
not fail the build. Run: python3 test/test_check_scenarios.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_scenarios as cs  # noqa: E402


def _write(d, calls_rows, truth_rows):
    calls = os.path.join(d, "calls.tsv")
    truth = os.path.join(d, "truth.tsv")
    out = os.path.join(d, "m.md")
    with open(calls, "w") as fh:
        fh.write("sample\tcaller\tsv_type\tbp5\tbp3\tsvlen\tsupport\thet\tcommon_deletion\textra\n")
        for r in calls_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(truth, "w") as fh:
        fh.write("#sample\tkind\tbp5\tbp3\tsvlen\thet\tdepth\texpect\n")
        for r in truth_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return calls, truth, out


def test_load_truth_reads_expect():
    d = tempfile.mkdtemp()
    _, truth, _ = _write(d, [], [("sv_del4977_h30", "del", 8469, 13447, 4977, 0.3, 300, "pass")])
    t = cs.load_truth(truth)
    assert t["sv_del4977_h30"][0] == ("del", 8469, 13447, "pass")


def test_main_renders_categorized_matrix_and_exits_zero():
    d = tempfile.mkdtemp()
    calls, truth, out = _write(
        d,
        # a wild-type FP (would gate if we gated — but this tool only reports)
        [("sv_wt", "eklipse", "deletion", 8470, 13447, 4977, 5, 0.01, 1, "x"),
         ("sv_del4977_h30", "mitohpc", "deletion", 8482, 13446, 4977, 99, 0.30, 1, "")],
        [("sv_del4977_h30", "del", 8469, 13447, 4977, 0.30, 300, "pass"),
         ("sv_dup", "dup", 6000, 7000, 0, 0.50, 300, "no_pass"),
         ("sv_inv_small", "inv", 6000, 6500, 501, 0.50, 300, "no_record"),
         ("sv_wt", "none", ".", ".", ".", 0, 300, "no_pass")])
    rc = cs.main(["--calls", calls, "--truth", truth, "--out-md", out,
                  "--samples", "sv_del4977_h30 sv_dup sv_inv_small sv_wt", "--no-real"])
    assert rc == 0
    md = open(out).read()
    assert "Accuracy (all scored scenarios)" in md
    assert "### Deletions" in md and "### Duplications" in md and "### Inversions" in md
    # eklipse made a common-del FP on wild-type -> reported (not gated)
    assert "FP" in md


def test_always_exits_zero_even_with_only_a_false_positive():
    d = tempfile.mkdtemp()
    calls, truth, out = _write(
        d, [("sv_wt", "eklipse", "deletion", 8470, 13447, 4977, 5, 0.01, 1, "x")],
        [("sv_wt", "none", ".", ".", ".", 0, 300, "no_pass")])
    rc = cs.main(["--calls", calls, "--truth", truth, "--out-md", out,
                  "--samples", "sv_wt", "--no-real"])
    assert rc == 0 and os.path.isfile(out)


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
