# Mitochondrial Structural Variation (mtDNA SV) caller suite

A batteries-included, end-to-end pipeline that runs **five mitochondrial
structural-variant / large-deletion callers** over a cohort of CRAM (or BAM)
files on an HPC cluster and consolidates their output into a single cohort
summary.

The callers are notoriously painful to install (Python 2, legacy samtools,
bundled binaries, dead reference-download URLs). This repo **Dockerizes all of
them into one image** built and published by GitHub Actions, and ships a simple
SLURM launcher that fans the work out with Apptainer.

| Caller | Method | Lang / runtime | Input it consumes here |
|--------|--------|----------------|------------------------|
| [eKLIPse](https://github.com/dooguypapua/eKLIPse) | soft-clip + BLAST breakpoints | Python 2.7, BLAST+, circos | realigned `chrM` BAM |
| [MitoSAlt](https://sourceforge.net/projects/mitosalt/) | LAST split-read clustering | Perl + R, LAST | mito FASTQ pair |
| [Splice-Break2](https://github.com/brookehjelm/Splice-Break2) | MapSplice2 junctions | bash + Python 2 + Java 8 | mito FASTQ pair |
| [MitoMut](https://github.com/shane-e945/MitoMut) | BLAT split-read | Python 3, pysam, BLAT | realigned `chrM` BAM |
| [MitoSeek](https://github.com/riverlee/MitoSeek) | discordant / large-TLEN reads | Perl, samtools 0.1.x | realigned `chrM` BAM |

> **Reference frame:** everything is normalised to the rCRS (NC_012920.1,
> contig `chrM`, 16569 bp). The classic ~4977 bp "common deletion"
> (m.8470_13447del) is the canonical positive control.

---

## How it works

```
            CRAM/BAM (any reference, any mito contig name)
                              │
              ┌───────────────▼───────────────┐
              │  preprocess.sh                 │
              │   detect mito contig (chrM/MT) │
              │   extract mito reads → FASTQ   │
              │   realign to rCRS (bwa mem)    │
              └───────┬───────────────┬────────┘
       chrM BAM (rCRS)│               │ mito FASTQ R1/R2
        ┌─────────────┼──────┐     ┌──┴───────────────┐
        ▼             ▼      ▼     ▼                  ▼
     eKLIPse       MitoMut  MitoSeek  MitoSAlt   Splice-Break2
        └─────────────┴──────┴─────┬──┴──────────────┘
                                   ▼
                         postprocess.py  →  cohort_*.tsv
```

* **One Docker image** ([`Dockerfile`](Dockerfile)) bundles all five callers,
  each isolated in its own micromamba environment (they have mutually
  incompatible dependencies — Python 2 vs 3, samtools 0.1.x vs 1.x, etc.).
* **One SLURM launcher** ([`slurm/run_mito_sv.sh`](slurm/run_mito_sv.sh)) takes a
  directory of CRAMs, fans out one array task per sample via Apptainer, then runs
  a dependent consolidation job.
* **All non-software prerequisites are handled automatically**: the rCRS
  reference + BLAST/LAST indexes are baked into the image; CRAM decoding uses a
  seeded reference cache (falling back to the EBI ENA MD5 service).

---

## Quick start (HPC)

```bash
# 1. Pull/point at the published image (or let the launcher pull it for you):
export MITO_SV_IMAGE=ghcr.io/jlanej/mitochondrial_structural_variation:latest

# 2. The only required argument is a directory of CRAMs (searched recursively):
./slurm/run_mito_sv.sh /path/to/crams

#   …optionally tune scheduler + outputs:
./slurm/run_mito_sv.sh /path/to/crams \
    --outdir   /scratch/$USER/mito_sv \
    --partition short --account my_alloc \
    --threads 8 --exts cram,bam
```

That submits:
1. a SLURM **array** (`mito-sv`) — one task per sample, each running all five
   callers inside the container;
2. a dependent **consolidation** job (`mito-sv-consolidate`) that writes the
   cohort summary once every sample finishes.

Results land under `--outdir`:

```
mito_sv_out/
├── samples.manifest.tsv
├── <sample>/
│   ├── preprocess/   <sample>.chrM.bam(.bai), <sample>.R1/R2.fastq.gz
│   ├── eklipse/ mitosalt/ splicebreak2/ mitomut/ mitoseek/
│   └── status.tsv
├── cohort_sv_calls.tsv          # every normalised call, long format
├── cohort_common_deletion.tsv   # sample × caller: common deletion detected?
├── cohort_caller_matrix.tsv     # sample × caller: deletion-call counts
└── cohort_summary.txt           # human-readable digest
```

### Run on a single node without a scheduler

```bash
./slurm/run_mito_sv.sh /path/to/crams --local      # loops samples here, no sbatch
./slurm/run_mito_sv.sh /path/to/crams --dry-run    # print what would run
```

---

## Cohort output schema

`cohort_sv_calls.tsv` (one row per call, per caller):

| column | meaning |
|--------|---------|
| `sample`, `caller` | identifiers |
| `sv_type` | `deletion` \| `duplication` \| `breakpoint` |
| `bp5`, `bp3` | 5′ / 3′ breakpoints (rCRS coordinates) |
| `svlen` | event size (bp) |
| `support` | supporting reads |
| `het` | heteroplasmy / allele fraction (0–1, normalised across callers) |
| `common_deletion` | 1 if it matches m.8470_13447del (±80 bp) |
| `extra` | caller-specific detail |

---

## Inputs & reference handling

* **Input:** CRAM or BAM, whole-genome or mito-only, any mito contig naming
  (`chrM` / `MT` / `M` / `NC_012920.1`). Detected automatically; override with
  `--mt-contig`.
* **CRAM decoding:** if you don't pass `--reference`, the pipeline uses a local
  `REF_CACHE` seeded with rCRS (so any rCRS-based CRAM decodes offline) and falls
  back to the EBI ENA MD5 service for other contigs (needs network). Pass
  `--reference genome.fa` to decode fully offline.
* **Normalisation:** mito reads are re-aligned to rCRS with `bwa mem`, giving
  every caller identical `chrM`/16569 input regardless of the source build (this
  also handles hg19's 16571 bp Yoruba mito and non-rCRS contig names).

---

## Building / CI

Two workflows:

* [`.github/workflows/docker-build.yml`](.github/workflows/docker-build.yml) —
  builds the image, runs the **functional smoke test against every caller** on
  the committed [test data](test/data), and (only if the smoke test passes, and
  only on `main`/tags) publishes to GHCR.
* [`.github/workflows/ci.yml`](.github/workflows/ci.yml) — fast checks on every
  push/PR: parser unit tests, shell + Python syntax, shellcheck.

The smoke test ([`test/smoke_test.sh`](test/smoke_test.sh)) mirrors the
[MitoHPC `sv-calling`](https://github.com/jlanej/MitoHPC/tree/sv-calling/test/sv)
test suite (the source of our test BAMs): it runs the full pipeline over the
whole diverse cohort and asserts —

* **Operating** — every caller runs to completion on the positive control and
  emits its expected output (hard gate);
* **Sensitivity** — the common deletion is detected by ≥1 caller on
  `sv_del4977_h30`, via the BAM **and** the CRAM input path (hard gate); plus a
  per-scenario sensitivity matrix across all constructs (common deletion at
  varying VAF/depth, non-repeat deletion, D-loop deletion, multi-deletion,
  duplication, origin-crossing, low coverage, real spike-in) — reported as
  warnings since the callers are heuristic;
* **Specificity** — no caller calls the common deletion on a sample that doesn't
  carry it (wild-type, duplication, origin-crossing, real healthy 1000G) — hard
  gate ([`test/check_scenarios.py`](test/check_scenarios.py));
* **Robustness** — degenerate inputs (wrong-contig / empty BAM) fail cleanly
  with no traceback (hard gate).

The scenario × caller matrix and per-caller operating/detection status are
written to [`test/example_output/SMOKE_SUMMARY.md`](test/example_output) on each
build. `SCOPE` controls breadth: `full` (default, all scenarios) on main,
`quick` (3-sample subset) on PRs.

Run the unit tests locally (no Docker needed):

```bash
python3 test/test_parsers.py            # caller output parsers
python3 test/test_check_scenarios.py    # scenario sensitivity/specificity gates
```

Run the full functional test against a locally built image:

```bash
docker build -t mito-sv:ci .          # heavy; CI normally does this
bash test/smoke_test.sh mito-sv:ci
```

> **Note:** the image is large (R, BLAST, circos, MapSplice2, four conda envs).
> Build it in CI / on a beefy node, not casually on a laptop.

### Test data

[`test/data`](test/data) vendors the MitoHPC `sv-calling` test corpus
(synthetic `chrM` BAMs with a [truth table](test/data/truth.tsv) plus a couple of
real 1000G chrM BAMs and a del4977 spike-in). The headline truth: del4977 lives
at breakpoints **8469–13447**.

---

## Repository layout

```
Dockerfile                 single image, 5 callers, 4 conda envs
docker/install/*.sh        per-caller install scripts (run during build)
assets/rCRS.chrM.fa        canonical rCRS reference (bundled)
vendor/MitoSAlt_1.1.1/     vendored MitoSAlt source (no upstream git repo)
pipeline/
  preprocess.sh            CRAM/BAM → normalised chrM BAM + mito FASTQ
  run_sample.sh            per-sample driver (preprocess + all callers)
  callers/run_*.sh         one wrapper per caller
  postprocess.py           cohort consolidation
  lib/parsers.py           per-caller output parsers (unit-tested)
slurm/
  run_mito_sv.sh           HPC launcher (the entry point)
  sample_job.sbatch        array task = one sample
  consolidate.sbatch       cohort summary job
test/
  test_parsers.py          parser unit tests
  smoke_test.sh            full functional CI against the image
  data/                    committed test BAMs + truth
```

---

## Limitations & notes

* **Linear rCRS alignment.** Realignment is to a linear rCRS; deletions spanning
  the artificial origin break (e.g. the `sv_origin` test case) may be missed by
  some callers. Most clinically relevant deletions (incl. the common deletion)
  are unaffected.
* **MitoSAlt runs in "enriched" mode** (LAST MT index only), which sidesteps the
  dead nuclear-genome download URLs in its `setup.sh`. NUMT discrimination is
  therefore lighter than a full nuclear+MT run.
* **Splice-Break2 / MitoSAlt realign from FASTQ**; eKLIPse / MitoMut / MitoSeek
  consume the normalised BAM. MitoSeek SV output for mito-only input is the
  large-deletion read set, which `parsers.py` clusters into calls.
* **Licensing.** eKLIPse is AGPL-3.0, MitoSeek GPL-2.0, MitoSAlt permissive;
  MitoMut and Splice-Break2 ship no explicit license and Splice-Break2 bundles
  MapSplice2 (academic-use). The published image is intended for research use —
  review each upstream license before redistribution.
* These callers are heuristic and were authored for targeted/enriched mtDNA
  data; treat cross-caller agreement (`cohort_caller_matrix.tsv`) as the signal,
  not any single caller in isolation.
