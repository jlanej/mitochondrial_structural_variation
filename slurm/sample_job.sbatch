#!/usr/bin/env bash
#SBATCH --job-name=mito-sv
###############################################################################
# sample_job.sbatch — one array task = one sample x ONE caller.
#
# Submitted by run_mito_sv.sh once per caller (MITO_SV_CALLER), --array=1-N with
# aftercorr on the prep array so each sample starts as soon as ITS prep is done.
# Consumes the prepared inputs (run_sample --prepared-bam) and writes to the
# isolated $OUT/by_caller/<caller>/<sample>/, so a slow caller (eKLIPse) never
# blocks a fast one and per-caller results land independently.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"
: "${MITO_SV_CALLER:?MITO_SV_CALLER (a single caller) required}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"; CALLER="$MITO_SV_CALLER"
STRICT="${MITO_SV_STRICT:-0}"
# Same packing as prep (identical NT) so caller chunk T covers the same samples
# as prep chunk T — which aftercorr:prep depends on, per chunk.
N="$(wc -l < "$MITO_SV_MANIFEST" | tr -d ' ')"
NT="${MITO_SV_NTASKS:-$N}"

if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: neither apptainer nor singularity on PATH" >&2; exit 1; fi

ndone=0; nfail=0
for ((L=TASK; L<=N; L+=NT)); do
    line="$(sed -n "${L}p" "$MITO_SV_MANIFEST")"
    [[ -n "$line" ]] || continue
    name="${line%%$'\t'*}"
    ndone=$((ndone+1))

    pdir="$OUT/prepared/$name"
    bam="$pdir/${name}.chrM.bam"; r1="$pdir/${name}.R1.fastq.gz"; r2="$pdir/${name}.R2.fastq.gz"
    if [[ ! -f "$bam" || ! -f "$r1" || ! -f "$r2" ]]; then
        echo "[$(date)] task=$TASK caller=$CALLER line=$L sample=$name: prepared inputs missing (prep failed?) — skip" >&2
        nfail=$((nfail+1)); continue
    fi

    sdir="$OUT/by_caller/$CALLER/$name"; mkdir -p "$sdir"
    cache="$OUT/.cache/${CALLER}_${L}"
    export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

    args=(--prepared-bam "$bam" --prepared-r1 "$r1" --prepared-r2 "$r2"
          --sample "$name" --outdir "$sdir" --threads "$THREADS" --callers "$CALLER")
    [[ "$STRICT" == "1" ]] && args+=(--strict)

    echo "[$(date)] task=$TASK caller=$CALLER line=$L sample=$name threads=$THREADS"
    if "$RT" exec --bind "$OUT" "$MITO_SV_SIF" /opt/pipeline/run_sample.sh "${args[@]}"; then
        echo "[$(date)] task=$TASK caller=$CALLER line=$L sample=$name finished"
    else
        echo "[$(date)] task=$TASK caller=$CALLER line=$L sample=$name FAILED" >&2; nfail=$((nfail+1))
    fi
done

# Non-strict (default): a caller failure never fails the task (matches per-sample
# behaviour — run_sample handles it internally). Strict: exit non-zero if any
# sample in the chunk failed, so the failure is visible.
echo "[$(date)] caller=$CALLER task=$TASK: $((ndone-nfail))/$ndone ok, $nfail failed"
if [[ "$STRICT" == "1" ]] && (( nfail > 0 )); then exit 1; fi
