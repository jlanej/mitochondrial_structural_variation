#!/usr/bin/env bash
#SBATCH --job-name=mito-extract
###############################################################################
# extract_job.sbatch — one array task = pull chrM out of ONE (whole-genome) CRAM,
# ONCE, idempotently.
#
# Real WGS CRAMs are large; this is the only stage that decodes them. The tiny
# chrM slice is cached at $OUT/chrM/<sample>.chrM.bam and reused by prep + every
# caller — and by re-runs: an already-valid slice is skipped, so the big CRAM is
# never re-touched. Runs as a parallel array up front; everything else depends on
# it (aftercorr), so a sample's downstream work starts as soon as ITS slice lands.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"; : "${MITO_SV_INPUT_DIR:?}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"
# Sample packing: this task handles manifest lines TASK, TASK+NT, TASK+2*NT, ...
# where NT is the array size (MITO_SV_NTASKS). NT defaults to N so a non-packed
# invocation (e.g. local one-per-task) still processes exactly one sample.
N="$(wc -l < "$MITO_SV_MANIFEST" | tr -d ' ')"
NT="${MITO_SV_NTASKS:-$N}"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

mkdir -p "$OUT/chrM"
ndone=0; nfail=0
for ((L=TASK; L<=N; L+=NT)); do
    line="$(sed -n "${L}p" "$MITO_SV_MANIFEST")"
    [[ -n "$line" ]] || continue
    name="${line%%$'\t'*}"; path="${line#*$'\t'}"
    ndone=$((ndone+1))
    out="$OUT/chrM/${name}.chrM.bam"

    # Per-sample micromamba cache (shared $HOME lock) — before any micromamba.
    cache="$OUT/.cache/extract_${L}"
    export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

    # Idempotent: a valid cached slice (bam + index, passes quickcheck) -> skip;
    # the CRAM is never re-decoded on a re-run.
    if [[ -f "$out" && -f "$out.bai" ]] && "$RT" exec --bind "$OUT" "$MITO_SV_SIF" \
            micromamba run -n mitosv samtools quickcheck "$out" 2>/dev/null; then
        echo "[$(date)] extract task=$TASK line=$L sample=$name: cached $out (skip)"
        continue
    fi

    work="$OUT/.cache/extract_work_${L}"; rm -rf "$work"; mkdir -p "$work"
    binds=(--bind "$MITO_SV_INPUT_DIR" --bind "$OUT")
    args=(--input "$path" --sample "$name" --outdir "$work" --threads "$THREADS" --extract-only "$out")
    if [[ -n "${MITO_SV_REFERENCE:-}" ]]; then
        binds+=(--bind "$(dirname "$MITO_SV_REFERENCE")"); args+=(--reference "$MITO_SV_REFERENCE")
    fi

    echo "[$(date)] extract task=$TASK line=$L sample=$name input=$path -> $out"
    if "$RT" exec "${binds[@]}" "$MITO_SV_SIF" /opt/pipeline/preprocess.sh "${args[@]}"; then
        echo "[$(date)] extract task=$TASK line=$L sample=$name done"
    else
        echo "[$(date)] extract task=$TASK line=$L sample=$name FAILED" >&2; nfail=$((nfail+1))
    fi
    rm -rf "$work"
done

# Isolated per-sample failures are logged but don't fail the task (so the chunk's
# other samples still flow downstream via aftercorr). Only a wholesale failure —
# every sample in the chunk — exits non-zero, surfacing a systemic problem.
echo "[$(date)] extract task=$TASK: $((ndone-nfail))/$ndone ok, $nfail failed"
if (( ndone > 0 && nfail == ndone )); then exit 1; fi
