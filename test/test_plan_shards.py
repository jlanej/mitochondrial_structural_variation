#!/usr/bin/env python3
"""Unit tests for the matrix shard planner (test/plan_shards.py)."""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plan_shards as P  # noqa: E402

TRUTH = "\n".join([
    "#sample\tkind\tbp5\tbp3\tsvlen\thet\tdepth\texpect",
    "sv_del4977_h30\tdel\t8469\t13447\t4977\t0.3\t300\tpass",
    "sv_del6000_h50\tdel\t5999\t10999\t4999\t0.5\t300\tpass",
    "sv_lowcov\tdel\t8469\t13447\t4977\t0.5\t40\tpass",
    "sv_dup\tdup\t6000\t7000\t0\t0.5\t300\tno_pass",
    "sv_inv_small\tinv\t6000\t6500\t501\t0.5\t300\tno_record",
    "sv_wt\tnone\t.\t.\t.\t0\t300\tno_pass",
]) + "\n"


def _truth_file():
    fd, p = tempfile.mkstemp(suffix=".tsv")
    os.write(fd, TRUTH.encode())
    os.close(fd)
    return p


def _all(shards):
    out = []
    for s in shards:
        out += s["samples"].split()
    return out


def test_every_sample_assigned_once():
    shards, _ = P.plan(_truth_file(), 3, "all")
    got = _all(shards)
    # 6 mock + 4 reals = 10 samples, each exactly once
    assert len(got) == 10 and len(set(got)) == 10
    assert "sv_wt" in got and "NA12718" in got and "spike_del4977_h20" in got


def test_extras_only_on_first_shard():
    shards, _ = P.plan(_truth_file(), 4, "all")
    extras = [s for s in shards if s["extras"]]
    assert len(extras) == 1 and extras[0]["idx"] == 0


def test_reals_are_spread_not_stacked():
    # the 4 deep real BAMs should land on distinct shards (LPT balances them)
    shards, _ = P.plan(_truth_file(), 4, "all")
    real_shards = [s["idx"] for s in shards
                   for x in s["samples"].split() if x in P.REALS]
    assert len(set(real_shards)) == 4


def test_balanced_loads():
    shards, loads = P.plan(_truth_file(), 4, "all")
    used = [loads[s["idx"]] for s in shards]
    # the heaviest shard is within ~1 real-weight of the lightest
    assert max(used) - min(used) <= P.REALS["NA12718"]


def test_del_suite_excludes_forward_looking():
    shards, _ = P.plan(_truth_file(), 3, "del")
    got = set(_all(shards))
    assert "sv_dup" not in got and "sv_inv_small" not in got
    assert "sv_del4977_h30" in got and "sv_wt" in got


def test_more_shards_than_samples_drops_empties():
    shards, _ = P.plan(_truth_file(), 50, "all")
    assert all(s["samples"].strip() for s in shards)   # no empty shards emitted


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
