#!/usr/bin/env python3
"""Generate one LOD grid cell: a chrM BAM carrying a deletion at a target
heteroplasmy (VAF) and depth, using MitoHPC's deterministic simulator + the
circular-aware realign (minimap2 -ax sr chrMC | -F 0x90C | circSam.pl | sort).

Mirrors lod_sweep.sim_run + realign verbatim (same mixing math, same per-cell
injective seed) so the reads are identical to MitoHPC's LOD harness. Output (in
--out): cell.r1.fastq.gz, cell.r2.fastq.gz, cell.bam(+.csi), truth.tsv.

Runs in the mitosv env (minimap2 + samtools + python3); circSam.pl needs perl.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import os
import random
import shutil
import subprocess
import sys
import tempfile
import zlib

RLEN, FMIN, FMAX, ERR, MTLEN = 150, 300, 450, 0.001, 16569


def load_simulator(path):
    spec = importlib.util.spec_from_file_location("make_testdata", path)
    mt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mt)
    return mt


def seed_for(variant, vaf, depth, rep):
    return zlib.crc32(("%s|%.4f|%d|%d" % (variant, vaf, depth, rep)).encode()) & 0x7FFFFFFF


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", required=True, help="deletion name, e.g. del4977")
    ap.add_argument("--bp5", type=int, required=True)
    ap.add_argument("--bp3", type=int, required=True)
    ap.add_argument("--vaf", type=float, required=True, help="target heteroplasmy fraction 0..1")
    ap.add_argument("--depth", type=int, required=True, help="outside-deletion per-base depth")
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--out", required=True, help="cell working dir")
    ap.add_argument("--refdir", default="/opt/MitoHPC/RefSeq",
                    help="dir with chrM.fa, chrMC.fa, chrM.fa.fai")
    ap.add_argument("--scriptsdir", default="/opt/MitoHPC/scripts", help="dir with circSam.pl")
    ap.add_argument("--simulator",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "make_testdata.py"))
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)

    mt = load_simulator(args.simulator)
    os.makedirs(args.out, exist_ok=True)
    ref = mt.read_fasta_single(os.path.join(args.refdir, "chrM.fa"))
    assert len(ref) == MTLEN, "chrM.fa is %d bp, expected %d" % (len(ref), MTLEN)

    seed = seed_for(args.variant, args.vaf, args.depth, args.rep)
    rng = random.Random(seed)

    # --- simulate the WT + deletion read mixture (lod_sweep.sim_run math) ---
    tmp = tempfile.mkdtemp(dir=args.out)
    r1, r2 = os.path.join(tmp, "r1.fq"), os.path.join(tmp, "r2.fq")
    with open(r1, "w") as f1, open(r2, "w") as f2:
        n_wt = round(args.depth * (1 - args.vaf) * len(ref) / RLEN)
        mt.emit_reads(f1, f2, ref + ref, n_wt, RLEN, FMIN, FMAX, ERR, rng, "wt")
        n_e = 0
        if args.vaf > 0:
            eg = mt.make_deletion(ref, args.bp5, args.bp3)
            n_e = round(args.depth * args.vaf * len(eg) / RLEN)
            mt.emit_reads(f1, f2, eg + eg, n_e, RLEN, FMIN, FMAX, ERR, rng, "ev")
    sys.stderr.write("[gen_cell] %s vaf=%.4f depth=%d rep=%d seed=%d : %d WT + %d event pairs\n"
                     % (args.variant, args.vaf, args.depth, args.rep, seed, n_wt, n_e))

    # --- circular-aware realign (MitoHPC's exact recipe) ---
    bam = os.path.join(args.out, "cell.bam")
    cmd = ("minimap2 -t %d -ax sr '%s/chrMC.fa' '%s' '%s' 2>/dev/null "
           "| samtools view -h -F 0x90C - "
           "| perl '%s/circSam.pl' -ref_len '%s/chrM.fa.fai' -offset 0 "
           "| samtools sort -@ %d -o '%s' --write-index - 2>/dev/null"
           % (args.threads, args.refdir, r1, r2, args.scriptsdir, args.refdir,
              args.threads, bam))
    rc = subprocess.run(["bash", "-c", cmd]).returncode
    if rc != 0 or not os.path.exists(bam):
        sys.stderr.write("[gen_cell] ERROR: realign failed (rc=%d)\n" % rc)
        return 1

    # --- gzip the FASTQ for the FASTQ-based callers, write truth ---
    for src, dst in ((r1, "cell.r1.fastq.gz"), (r2, "cell.r2.fastq.gz")):
        with open(src, "rb") as i, gzip.open(os.path.join(args.out, dst), "wb") as o:
            shutil.copyfileobj(i, o)
    shutil.rmtree(tmp, ignore_errors=True)

    with open(os.path.join(args.out, "truth.tsv"), "w") as fh:
        fh.write("#variant\tbp5\tbp3\tsvlen\tvaf\tdepth\trep\tseed\n")
        fh.write("%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\n"
                 % (args.variant, args.bp5, args.bp3, args.bp3 - args.bp5 - 1,
                    args.vaf, args.depth, args.rep, seed))
    sys.stderr.write("[gen_cell] wrote %s\n" % bam)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
