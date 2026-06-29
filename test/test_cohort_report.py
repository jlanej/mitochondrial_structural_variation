#!/usr/bin/env python3
"""Unit tests for the cohort summary report (pipeline/cohort_report.py):
sample discovery across both layouts (no double-count), the call filter, and the
PASS derivation."""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lib"))
import cohort_report as CR  # noqa: E402


def _mixed_tree():
    """sampleX in BOTH classic and by_caller; sampleY classic-only; sampleZ
    by_caller-only; plus noise dirs that must be ignored."""
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "sampleX", "mitohpc"))
    open(os.path.join(d, "sampleX", "status.tsv"), "w").close()
    os.makedirs(os.path.join(d, "by_caller", "eklipse", "sampleX", "eklipse"))
    os.makedirs(os.path.join(d, "sampleY", "mitohpc"))
    open(os.path.join(d, "sampleY", "status.tsv"), "w").close()
    os.makedirs(os.path.join(d, "by_caller", "mitosalt", "sampleZ", "mitosalt"))
    os.makedirs(os.path.join(d, "prepared", "sampleX"))
    os.makedirs(os.path.join(d, "logs"))
    return d


def test_no_double_count_mixed_layout():
    pairs, _ = CR.discover(_mixed_tree())
    names = [n for n, _ in pairs]
    assert names.count("sampleX") == 1, "sampleX double-counted across layouts"
    # the retained sampleX pair is the by_caller one (current run wins)
    sx = [sd for n, sd in pairs if n == "sampleX"][0]
    assert "by_caller" in sx


def test_keeps_both_classic_only_and_bycaller_only():
    pairs, _ = CR.discover(_mixed_tree())
    names = {n for n, _ in pairs}
    assert "sampleY" in names and "sampleZ" in names


def test_noise_dirs_excluded():
    pairs, _ = CR.discover(_mixed_tree())
    names = {n for n, _ in pairs}
    assert "prepared" not in names and "logs" not in names and "by_caller" not in names


def test_present_callers_detected():
    _, present = CR.discover(_mixed_tree())
    # eklipse + mitosalt have by_caller dirs; mitohpc has a classic subdir
    assert {"eklipse", "mitosalt", "mitohpc"} <= present


def test_is_call_filters_raw_evidence():
    # a typed deletion with both breakpoints is a call
    assert CR._is_call({"sv_type": "deletion", "bp5": 8470, "bp3": 13447})
    # MitoSeek's raw discordant 'breakpoint' (no bp3) is NOT a call
    assert not CR._is_call({"sv_type": "breakpoint", "bp5": 8470, "bp3": None})
    # a deletion missing a breakpoint is not a call
    assert not CR._is_call({"sv_type": "deletion", "bp5": 8470, "bp3": None})


def test_passed_derivation():
    # MitoHPC: PASS only when its FILTER is PASS
    assert CR._passed({"caller": "mitohpc", "extra": "filter=PASS;flags=."})
    assert not CR._passed({"caller": "mitohpc", "extra": "filter=lowJR;flags=."})
    # other callers: their output IS the final call set -> always pass
    assert CR._passed({"caller": "eklipse", "extra": ""})


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
