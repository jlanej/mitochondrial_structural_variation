#!/usr/bin/env bash
#SBATCH --job-name=mito-prep
###############################################################################
# prep_job.sbatch — one array task = realign ONE sample's cached chrM slice to
# rCRS -> normalised chrM BAM + mito FASTQ in $OUT/prepared/<sample>/.
#
# Reads the small $OUT/chrM/<sample>.chrM.bam produced by the mito-extract stage
# (NOT the original WGS CRAM — that is decoded only once, in extract). The
# per-caller arrays then consume these prepared inputs via run_sample
# --prepared-bam, so the CRAM is never re-touched.
###############################################################################
set -euo pipefail
: "${MITO_SV_SIF:?}"; : "${MITO_SV_OUTDIR:?}"; : "${MITO_SV_MANIFEST:?}"
THREADS="${MITO_SV_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
TASK="${SLURM_ARRAY_TASK_ID:?must run as a SLURM array job}"
OUT="$MITO_SV_OUTDIR"
# Same packing as extract (identical NT) so prep chunk T covers exactly the same
# samples as extract chunk T — which aftercorr:extract depends on, per chunk.
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
    chrm="$OUT/chrM/${name}.chrM.bam"
    if [[ ! -f "$chrm" ]]; then
        echo "[$(date)] prep task=$TASK line=$L sample=$name: chrM missing ($chrm; extract failed?) — skip" >&2
        nfail=$((nfail+1)); continue
    fi
    pdir="$OUT/prepared/$name"; mkdir -p "$pdir"

    # Per-sample micromamba cache: Apptainer binds the shared host $HOME, so without
    # isolation every concurrent micromamba serialises on one mamba/proc lock.
    cache="$OUT/.cache/prep_${L}"
    export APPTAINERENV_XDG_CACHE_HOME="$cache" SINGULARITYENV_XDG_CACHE_HOME="$cache"

    echo "[$(date)] prep task=$TASK line=$L sample=$name input=$chrm threads=$THREADS"
    if "$RT" exec --bind "$OUT" "$MITO_SV_SIF" /opt/pipeline/preprocess.sh \
            --input "$chrm" --sample "$name" --outdir "$pdir" --threads "$THREADS"; then
        echo "[$(date)] prep task=$TASK line=$L sample=$name -> $pdir"
    else
        echo "[$(date)] prep task=$TASK line=$L sample=$name FAILED" >&2; nfail=$((nfail+1))
    fi
done

echo "[$(date)] prep task=$TASK: $((ndone-nfail))/$ndone ok, $nfail failed"
if (( ndone > 0 && nfail == ndone )); then exit 1; fi
