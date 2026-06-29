#!/usr/bin/env python3
"""Unit tests for the MitoBreak known-breakpoint matcher (pipeline/lib/mitobreak.py)
against the vendored assets/mitobreak.tsv.gz."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lib"))
import mitobreak as MB  # noqa: E402

DB = MB.load()


def test_db_loads_with_both_types():
    assert DB, "MitoBreak DB failed to load (assets/mitobreak.tsv.gz missing?)"
    assert len(DB.get("DEL", [])) > 1000, "expected >1000 catalogued deletions"
    assert len(DB.get("DUP", [])) > 10, "expected catalogued duplications"


def test_common_deletion_is_known():
    # the ~4977 bp common deletion (our convention bp5=8470, bp3=13447)
    m = MB.match(DB, "deletion", 8470, 13447)
    assert m is not None, "common deletion should match a MitoBreak entry"
    assert m[0].startswith("DEL_"), "match id should be a deletion id"


def test_off_target_call_is_novel():
    # a deletion nowhere near any catalogued breakpoint -> novel
    assert MB.match(DB, "deletion", 200, 400, tol=20) is None


def test_tolerance_is_enforced():
    # 8470/13447 matches; shifting both breakpoints well beyond tol must not match
    assert MB.match(DB, "deletion", 8470 + 100, 13447 + 100, tol=20) is None


def test_inversion_never_matches():
    # the catalogue has no inversions; an inversion call is always novel
    assert MB.match(DB, "inversion", 8470, 13447) is None


def test_missing_breakpoints_return_none():
    assert MB.match(DB, "deletion", None, 13447) is None
    assert MB.match(DB, "deletion", 8470, None) is None


def test_closest_entry_wins():
    # exact-ish hit returns the nearest catalogued id (summed distance small)
    m = MB.match(DB, "deletion", 8470, 13447)
    assert m[3] <= 4, "nearest common-deletion entry should be within a few bp"


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
