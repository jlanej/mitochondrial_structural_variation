#!/usr/bin/env python3
"""Plan balanced scenario shards for the matrix CI.

Each shard runs a subset of samples (all callers) on its own runner. Cost is
dominated by the slow callers (eKLIPse / Splice-Break2), which scale with read
count, i.e. depth — so we weight each sample by its depth (the real 1000G BAMs
at ~2400x are ~8x a 300x mock) and assign greedily by longest-processing-time
(LPT): heaviest sample to the currently-lightest shard. That spreads the
expensive deep samples across shards and keeps the slow callers off any single
shard's critical path.

Prints a JSON array for the GitHub Actions matrix:
  [{"idx":0,"samples":"sv_a sv_b ...","extras":"1"}, ...]
Shard 0 carries the cheap "extras" (CRAM round-trip + degenerate-input checks).
"""
import argparse
import json

# Committed real 1000G BAMs + the spike — high depth (not in truth.tsv).
REALS = {"spike_del4977_h20": 2400.0, "NA12718": 2400.0,
         "NA12748": 2400.0, "NA12775": 2400.0}
DEFAULT_DEPTH = 300.0


def load_truth(path):
    depth, kinds = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            s = f[0]
            kinds.setdefault(s, set()).add(f[1])
            try:
                d = float(f[6])
            except (IndexError, ValueError):
                d = DEFAULT_DEPTH
            # multi-event samples (e.g. multidel) cost a bit more per depth-unit.
            depth[s] = depth.get(s, 0.0) + max(d, 1.0)
    return depth, kinds


def plan(truth, n_shards, suite="all"):
    depth, kinds = load_truth(truth)
    samples = []
    for s in sorted(kinds):
        if suite == "del" and not (kinds[s] & {"del", "delwrap", "none"}):
            continue
        samples.append(s)
    samples += list(REALS)
    weight = {s: REALS.get(s, depth.get(s, DEFAULT_DEPTH)) for s in samples}

    n = max(1, min(n_shards, len(samples)))
    shards = [[] for _ in range(n)]
    loads = [0.0] * n
    for s in sorted(samples, key=lambda x: -weight[x]):
        i = loads.index(min(loads))
        shards[i].append(s)
        loads[i] += weight[s]

    out = []
    for i, sh in enumerate(shards):
        if not sh:
            continue
        out.append({"idx": i, "samples": " ".join(sh),
                    "extras": "1" if i == 0 else ""})
    return out, loads


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--shards", type=int, default=6)
    ap.add_argument("--suite", default="all", choices=["all", "del"])
    args = ap.parse_args(argv)
    shards, loads = plan(args.truth, args.shards, args.suite)
    print(json.dumps(shards))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
