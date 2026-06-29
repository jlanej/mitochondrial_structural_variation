# Mitochondrial Structural Variation (mtDNA SV) caller suite

A batteries-included, end-to-end pipeline that runs **six mitochondrial
structural-variant / large-deletion callers** over a cohort of CRAM (or BAM)
files on an HPC cluster, consolidates their output into a single cohort summary,
and produces an **interactive caller-comparison report**.

The callers are notoriously painful to install (Python 2, legacy samtools,
bundled binaries, dead reference-download URLs). This repo **Dockerizes all of
them into one image** built and published by GitHub Actions, and ships a simple
SLURM launcher that fans the work out with Apptainer.

| Caller | Method | Lang / runtime | Input it consumes here |
|--------|--------|----------------|------------------------|
| [MitoHPC](https://github.com/jlanej/MitoHPC/tree/sv-calling) *(reference)* | split-read + coverage-drop | Python 3, pysam | realigned `chrM` BAM |
| [eKLIPse](https://github.com/dooguypapua/eKLIPse) | soft-clip + BLAST breakpoints | Python 2.7, BLAST+, circos | realigned `chrM` BAM |
| [MitoSAlt](https://sourceforge.net/projects/mitosalt/) | LAST split-read clustering | Perl + R, LAST | mito FASTQ pair |
| [Splice-Break2](https://github.com/brookehjelm/Splice-Break2) | MapSplice2 junctions | bash + Python 2 + Java 8 | mito FASTQ pair |
| [MitoMut](https://github.com/shane-e945/MitoMut) | BLAT split-read | Python 3, pysam, BLAT | realigned `chrM` BAM |
| [MitoSeek](https://github.com/riverlee/MitoSeek) | discordant / large-TLEN reads | Perl, samtools 0.1.x | realigned `chrM` BAM |

[**MitoHPC**](https://github.com/jlanej/MitoHPC/tree/sv-calling)'s own SV caller
is bundled as the **reference method**: it runs on the same normalised `chrM`
BAM as the others, giving an apples-to-apples comparison. See the live
[comparison report](docs/index.html) (`docs/index.html`) — caller × scenario
detection matrix, runtime, and sensitivity, regenerated on every build.

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
  directory of CRAMs, preprocesses each once (`mito-prep`), then fans out **one
  array per caller** (so a slow caller never blocks a fast one) and refreshes the
  cohort summary as each caller finishes — see below.
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

That submits a **per-caller** pipeline so no single slow caller (e.g. eKLIPse) can
hold up the rest, and results stream in as each caller completes. The compute jobs
(extract / prep / callers) use one generous profile — **24 threads, 64 GB, 10 h
walltime** (tune with `--threads / --mem / --time`); the lightweight consolidations
run on **2 threads / 16 GB / 2 h** (`MITO_SV_CONS_{CPUS,MEM,TIME}`):
1. **`mito-extract`** — pull chrM out of the (whole-genome) CRAM, **once**, into
   `$OUT/chrM/<sample>.chrM.bam`. This is the only stage that decodes the big CRAM,
   and it is **idempotent** — an already-valid slice is skipped, so re-runs never
   re-touch the CRAM.
2. **`mito-prep`** — `aftercorr` on extract: realign each cached slice to rCRS → chrM
   BAM + FASTQ in `prepared/<sample>/` (works off the small BAM, not the CRAM).
3. **`mito-<caller>`** — one array *per caller*, `aftercorr` on `mito-prep` (a chunk's
   caller starts the moment *its* prep is done). Each writes to an isolated
   `by_caller/<caller>/<sample>/`.

Every stage is submitted as a **fixed number of array tasks** (`--chunks`, default
100) rather than one task per sample: each task *strides* over the manifest and
processes its slice of the cohort sequentially (`ceil(n/chunks)` samples each). So
the submitted-job count is **constant** (`~8·chunks + 7`, e.g. ~807 for six callers)
no matter how large the cohort grows — a 2 000-sample run submits 807 jobs, not
16 007, staying well under `AssocMaxSubmitJobLimit`. The trade-off is walltime: a
task runs its samples back-to-back, so raise `--time` for big cohorts (or raise
`--chunks` to keep slices small, as long as `8·chunks + 7` stays under the limit —
the launcher refuses to submit otherwise). A single sample failing inside a chunk is
logged and skipped; the chunk's other samples still flow downstream.
4. **`cons-<caller>`** — fires after each caller finishes **all** its CRAMs:
   rebuilds the cohort tables **and** the interactive
   [`cohort_sv_summary.html`](docs/cohort_sv_summary.html), so the summary appears
   as soon as the first caller completes and fills in as the rest land.
5. **`mito-final`** — after everything: the authoritative summary, then drops the
   bulky `prepared/` inputs (the `chrM/` slices are kept for fast idempotent re-runs).

Concurrent cells/consolidations each get an isolated `XDG_CACHE_HOME` (Apptainer
binds the shared host `$HOME`, so otherwise every `micromamba run` serialises on
one cache lock); consolidations are flock-serialised. Running tasks stay bounded by
a global concurrency budget (`--max-concurrent`, default 500) that is split evenly
across the caller arrays — so their combined `%N` throttles sum to ~500 rather than
500 per caller, which keeps the cluster busy without tripping the scheduler's
submit-job limit.

Results land under `--outdir`:

```
mito_sv_out/
├── samples.manifest.tsv
├── chrM/<sample>.chrM.bam(.bai)  # chrM sliced from the CRAM once (kept; idempotent)
├── prepared/<sample>/            # realigned chrM BAM + FASTQ (removed by mito-final)
├── by_caller/<caller>/<sample>/  # <caller>/ outputs + status.tsv
├── cohort_sv_calls.tsv           # every normalised call, long format
├── cohort_common_deletion.tsv    # sample × caller: common deletion detected?
├── cohort_caller_matrix.tsv      # sample × caller: deletion-call counts
├── cohort_summary.txt            # human-readable digest
└── cohort_sv_summary.html        # interactive global-metrics report (see below)
```

### Cohort summary report (`cohort_sv_summary.html`)

A self-contained, interactive report of **global** caller-comparison metrics for a
grant reviewer — no per-call table, no samplot, nothing per-individual that would
bloat it. Generated by [`pipeline/cohort_report.py`](pipeline/cohort_report.py):

* **known vs novel per caller** — the % of each caller's PASS calls whose breakpoint
  is catalogued in **MitoBreak** (vendored [`assets/mitobreak.tsv.gz`](assets);
  matched within 20 bp on both ends, the same rule MitoHPC uses), a scale-free
  credibility/specificity proxy;
* **SVs called per individual** — a box plot per caller (log axis — callers differ
  by orders of magnitude), all calls vs PASS-only;
* **calls per individual by SV type** — type composition + a mean-per-individual
  table;
* **common deletion (del4977) detection** and **reported-heteroplasmy** distributions.

A *call* is a typed SV (deletion/duplication/inversion) with both breakpoints; raw
discordant-read evidence is not counted.

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
  builds + publishes the image and runs the **functional smoke test against every
  caller** on the committed [test data](test/data). On `main`/tags the scenario
  cohort is **fanned out across a matrix** of runner jobs: a planner
  ([`test/plan_shards.py`](test/plan_shards.py)) splits the samples into balanced
  shards (LPT bin-packing weighted by depth, so the deep real 1000G BAMs spread
  out and the slow callers stay off any single shard's critical path), each shard
  runs its subset serially on its own runner (no in-process concurrency → the
  memory-hungry callers don't OOM each other), and a **consolidate** job merges
  the shard outputs, post-processes, gates, and builds the report. PRs run a quick
  single-job smoke. `smoke_test.sh` is phase-aware (`full` / `shard` /
  `consolidate`) so all three reuse the same logic.
* [`.github/workflows/ci.yml`](.github/workflows/ci.yml) — fast checks on every
  push/PR: parser unit tests, shell + Python syntax, shellcheck.

The smoke test ([`test/smoke_test.sh`](test/smoke_test.sh)) runs the full
pipeline over the whole [MitoHPC `sv-calling`](https://github.com/jlanej/MitoHPC/tree/sv-calling/test/sv)
test cohort (the source of our test BAMs). It separates two concerns:

* **Hard gates** — only things the *pipeline* controls:
  - **Operating** — every caller runs to completion on the positive control
    `sv_del4977_h30` and emits its expected output;
  - post-processing produces the cohort tables;
  - **Robustness** — degenerate inputs (wrong-contig / empty BAM) fail cleanly,
    no traceback.
* **Evaluation only (never gates)** — how each *third-party caller* behaves
  across the diverse constructs (common deletion at varying VAF/depth, non-repeat
  deletion, D-loop deletion, multi-deletion, duplication, origin-crossing, low
  coverage, real spike-in / healthy 1000G). We don't control the callers' source,
  so their sensitivity/specificity is **recorded as a caller-comparison matrix**
  ([`test/check_scenarios.py`](test/check_scenarios.py)), not asserted.

The scenario × caller comparison matrix and per-caller operating/detection
status are written to
[`test/example_output/SMOKE_SUMMARY.md`](test/example_output) on each build.
`SCOPE` controls breadth: `full` (default, all scenarios) on main, `quick`
(3-sample subset) on PRs.

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

[`test/data`](test/data) vendors the MitoHPC `sv-calling` test corpus — **21
synthetic `chrM` BAMs** spanning deletions (incl. the common del4977 at
**8469–13447**, a non-repeat del, D-loop, multi-deletion, near-homoplasmy,
low-coverage, and the 45 bp / 500 bp / 13 kb size bounds), origin-crossing
deletions, tandem duplications, inversions, and complex (dup-del / inverted-dup)
events — plus three real 1000G chrM BAMs and a del4977 spike-in. The per-sample
`kind`/`expect` labels live in [truth.tsv](test/data/truth.tsv) and the
categories + hover descriptions in [scenarios.json](test/data/scenarios.json).
Caller behaviour on every scenario is **evaluation only** (the deletion/dup/inv
expectations are MitoHPC's assertions about *its* caller, not gates for ours).

The CI scenario suite runs the full set by default; flip `SUITE=del` in
[docker-build.yml](.github/workflows/docker-build.yml) (or pass `del` as the 4th
arg to `test/smoke_test.sh`) to run only the **deletions + controls** (13 BAMs,
skipping the forward-looking dup/inv/complex BAMs) if the full suite is too slow.

The interactive report scores each caller as a binary deletion detector and
breaks results down **by category** — the detection matrix defaults to the
**Deletions** tab with tabs for each category plus **All**, and reports
per-category sensitivity, specificity, precision, F1, balanced accuracy and MCC
(shared kernel: [pipeline/lib/sv_eval.py](pipeline/lib/sv_eval.py), unit-tested).
Hover any BAM name in the matrix for what that scenario tests.

---

## Limit-of-detection (LOD) sweep

A second batteries-included HPC entry point benchmarks **how low a heteroplasmy
each caller can detect**, modeled on MitoHPC's own LOD methods but run for all
six callers in parallel.

```bash
./slurm/run_lod.sh                     # full grid (default)
./slurm/run_lod.sh --quick             # small grid
./slurm/run_lod.sh --hets 0,0.05,0.1 --depths 1000,2000 --reps 10 --deletions del4977
```

For every cell `(deletion ∈ {del4977, del6000} × heteroplasmy × depth ×
replicate)` it simulates a chrM BAM carrying that deletion at that VAF/depth
(MitoHPC's deterministic [`make_testdata.py`](pipeline/lod/make_testdata.py),
injective per-cell seed), runs all callers under **two input arms** — `pipeline`
(bwa-mem normalization, production behaviour) and `circular` (MitoHPC's
`minimap2 + circSam.pl` circular-aware BAM) — and scores detection (a call within
**30 bp** summed breakpoint error). The default grid is **10 heteroplasmies × 4
depths (250–2000×) × 10 replicates × 2 deletions = 800 cells**; heteroplasmy is
dense at the low end (where the LOD lives) and depth tops out at a real-world
2000× to expose how depth warps runtime.

Work is split **per (caller, depth)**: one SLURM array per combination, named
`lod-<caller>-d<depth>`, each running only its caller so no single slow caller
(eKLIPse) can monopolise a submission. Within a task, cells run concurrently
across the 24 cores (`--cpus-per-task 24`, `THREADS/TPC = 24/2 = 12` at once); the
800-cell grid fans out to ~408 tasks. Jobs are uncapped unless the projected total
exceeds `--max-jobs` (500), at which point cells-per-job is raised (iterations
distributed) until it fits, one caller per job. **Consolidation cascades** as
results land — `cons-<caller>-d<depth>` after each array (prints that scope's
runtime + a rough remaining-wall projection, the realtime cancel signal),
`cons-<caller>` after a caller's depths, and `lod-final` after everything — so a
cumulative report appears as soon as the first caller finishes a depth.

Each concurrent cell gets its own `XDG_CACHE_HOME` (Apptainer binds the shared
host `$HOME`, so otherwise every `micromamba run` cluster-wide serializes on one
`mamba/proc` lock — the failure mode that once stalled the whole sweep). A
background **heartbeat** logs each task's `done/total` plus every active cell's
elapsed time, phase, and finished callers (so a hung cell/caller is obvious), and
heavy intermediates are deleted the moment a cell is scored. The consolidations
aggregate the sweep:

* per `(caller, depth, deletion, arm)`: a detection-probability curve over
  heteroplasmy → **LOD50 / LOD95** (Firth-penalized logistic, separation-robust;
  pure-Python, no scipy) with Wilson + cluster-bootstrap CIs and an empirical
  transition read-out;
* a **runtime summary** (`lod_runtime.tsv`) — per-cell wall-clock distribution
  **broken down by depth** (and per arm): n, mean, median, p25/p75, min/max,
  total seconds, so the depth→cost scaling is explicit;
* an interactive **`lod_report/index.html`** — methods, per-caller LOD curves,
  detection heatmap, an LOD summary table, **runtime box-and-whisker plots + a
  runtime-vs-depth chart and per-depth table**, pipeline-vs-circular comparison,
  and an interpretation guide.

All callers (including the MitoHPC reference) run in the **same image, env, and
node**, each timed identically by `run_sample.sh` (`t0=$SECONDS … $((SECONDS-t0))`),
so the runtimes are a fair head-to-head; the `.sif` is built once per submission
and reused across array tasks (no per-run image pull). MitoHPC's optional samplot
visualization is explicitly disabled (`HP_SV_PLOT=`) so its runtime is pure SV
calling. Because MitoHPC is our own repo, the image **tracks its `sv-calling`
branch**: CI resolves the latest commit and passes it as the `MITOHPC_REF` build
arg, so changes are benchmarked on the next build (the exact commit is recorded at
`/opt/MitoHPC/GIT_SHA`). The third-party callers stay pinned for reproducibility.

CI runs a **single-iteration gate** ([`test/lod_smoke.sh`](test/lod_smoke.sh)) —
one tiny cell through both arms — to prove the machinery runs; the full sweep is
the HPC job. The LOD statistics are unit-tested
([`test/test_lod_stats.py`](test/test_lod_stats.py), no Docker).

---

## Repository layout

```
Dockerfile                 single image, 6 callers, 4 conda envs
docker/install/*.sh        per-caller install scripts (run during build)
assets/rCRS.chrM.fa        canonical rCRS reference (bundled)
vendor/MitoSAlt_1.1.1/     vendored MitoSAlt source (no upstream git repo)
pipeline/
  preprocess.sh            CRAM/BAM → normalised chrM BAM + mito FASTQ
  run_sample.sh            per-sample driver (preprocess + all callers; --prepared)
  callers/run_*.sh         one wrapper per caller (incl. run_mitohpc.sh)
  postprocess.py           cohort consolidation (+ cohort_runtime.tsv)
  make_report.py           interactive docs/index.html generator (scenario suite)
  cohort_report.py         interactive cohort_sv_summary.html (real-cohort metrics)
  lib/parsers.py           per-caller output parsers (unit-tested)
  lib/sv_eval.py           categorize scenarios + per-category accuracy metrics
  lib/mitobreak.py         MitoBreak known-breakpoint matcher (unit-tested)
  lod/                     limit-of-detection sweep tooling
    make_testdata.py       vendored MitoHPC read simulator
    gen_cell.py            simulate one (deletion,vaf,depth) cell -> chrM BAM
    run_cell.sh            one cell: generate + run callers (both arms) + score
    score_cell.py          score detection vs truth (BP_TOL=30)
    lod_stats.py           pure-Python LOD stats kernel (Firth logistic, Wilson)
    analyze_lod.py         sweep -> lod_cells.tsv + lod_fits.tsv + lod_runtime.tsv
    make_lod_report.py     interactive LOD report generator
slurm/
  run_mito_sv.sh           cohort HPC launcher (extract-once -> per-caller arrays)
  extract_job.sbatch       array task = slice chrM from one CRAM, once, idempotent
  prep_job.sbatch          array task = realign one sample's cached chrM slice to rCRS
  sample_job.sbatch        array task = one sample x ONE caller (prepared inputs)
  consolidate.sbatch       cohort tables + cohort_sv_summary.html; per-caller + final
  run_lod.sh               LOD-sweep launcher (per-(caller,depth) arrays + cascade)
  lod_array.sbatch         LOD array task = a chunk of cells for ONE caller (24-thread)
  lod_consolidate.sbatch   LOD analyze + report; fires per-(caller,depth), per-caller, final
test/
  test_parsers.py          parser unit tests
  test_sv_eval.py          scenario categorization + accuracy-metric unit tests
  test_check_scenarios.py  scenario-evaluator unit tests
  test_lod_stats.py        LOD statistics unit tests
  check_scenarios.py       truth-driven caller-comparison evaluator
  smoke_test.sh            functional CI; phase-aware full | shard | consolidate
  plan_shards.py           balance scenario samples into matrix shards (LPT)
  lod_smoke.sh             single-iteration LOD CI gate
  data/                    committed test BAMs + truth.tsv + scenarios.json
docs/index.html            interactive caller-comparison report (CI-generated)
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
* **MitoHPC calls are always typed `deletion`.** `parse_mitohpc` (in `parsers.py`)
  labels every row `sv_type=deletion` and does not read MitoHPC's own `svtype`
  column. Harmless as run here — MitoHPC runs in its default deletion-only mode, so
  every row already is a deletion — but if it is ever run with `HP_SV_DUP` /
  `HP_SV_INV` enabled, its `DUP` / `INV` calls would be mislabeled `deletion` in
  `cohort_sv_calls.tsv`. Binary-detector scoring is unaffected (a span call on a
  duplication already counts as a deletion-like call); supporting those event types
  would mean mapping the `svtype` column in `parse_mitohpc`.
* **Licensing.** eKLIPse is AGPL-3.0, MitoSeek GPL-2.0, MitoSAlt permissive;
  MitoMut and Splice-Break2 ship no explicit license and Splice-Break2 bundles
  MapSplice2 (academic-use). The published image is intended for research use —
  review each upstream license before redistribution.
* These callers are heuristic and were authored for targeted/enriched mtDNA
  data; treat cross-caller agreement (`cohort_caller_matrix.tsv`) as the signal,
  not any single caller in isolation.
