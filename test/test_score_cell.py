#!/usr/bin/env python3
"""Unit tests for score_cell.py row emission.

Guards the per-caller LOD split: score_cell must emit one row per caller that
ACTUALLY RAN (per status.tsv), never a spurious detected=0 row for callers that
were not invoked. Before the fix it iterated all six callers unconditionally, so a
single-caller arm produced 5 bogus zero rows that collapsed detection rates.
"""
import csv
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lod"))
import score_cell as SC  # noqa: E402

TRUTH = ("variant\tbp5\tbp3\tsvlen\tvaf\tdepth\trep\tseed\n"
         "del4977\t8469\t13447\t4977\t0.05\t2000\t0\t0\n")


def _run(status_rows):
    """Score an empty sample dir whose status.tsv lists `status_rows`
    [(caller, status, seconds), ...]; return the emitted data rows (dicts)."""
    d = tempfile.mkdtemp()
    sd = os.path.join(d, "sample"); os.makedirs(sd)
    with open(os.path.join(sd, "status.tsv"), "w") as fh:
        fh.write("caller\tstatus\tseconds\n")
        for c, st, secs in status_rows:
            fh.write("%s\t%s\t%s\n" % (c, st, secs))
    truth = os.path.join(d, "truth.tsv")
    with open(truth, "w") as fh:
        fh.write(TRUTH)
    out = os.path.join(d, "shard.tsv")
    SC.main(["--sample-dir", sd, "--truth", truth, "--arm", "pipeline",
             "--sample", "s", "--out", out])
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_single_caller_emits_one_row():
    rows = _run([("eklipse", "ok", "46")])
    assert len(rows) == 1, "expected 1 row, got %d" % len(rows)
    assert rows[0]["caller"] == "eklipse"
    assert rows[0]["status"] == "ok"


def test_all_callers_emit_six_rows():
    rows = _run([(c, "ok", "5") for c in
                 ["mitohpc", "eklipse", "mitosalt", "splicebreak2", "mitomut", "mitoseek"]])
    assert len(rows) == 6, "expected 6 rows, got %d" % len(rows)
    assert {r["caller"] for r in rows} == set(SC.parsers.CALLERS)


def test_caller_that_ran_but_failed_still_scored():
    # a caller that ran and crashed is 'not detected' (detected=0), not omitted
    rows = _run([("mitosalt", "failed", "3")])
    assert len(rows) == 1 and rows[0]["caller"] == "mitosalt"
    assert rows[0]["detected"] == "0" and rows[0]["status"] == "failed"


def test_no_callers_ran_emits_no_data_rows():
    rows = _run([])
    assert len(rows) == 0, "expected 0 rows, got %d" % len(rows)


def test_row_order_follows_caller_catalogue():
    rows = _run([("mitoseek", "ok", "5"), ("mitohpc", "ok", "5")])
    # emitted in parsers.CALLERS order regardless of status.tsv order
    assert [r["caller"] for r in rows] == ["mitohpc", "mitoseek"]


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
