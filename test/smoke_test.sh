#!/usr/bin/env bash
###############################################################################
# smoke_test.sh — robust functional CI for the built image.
#
#   bash test/smoke_test.sh [IMAGE] [CALLERS]
#
# Runs the full per-sample pipeline (preprocess + all callers) inside the image
# on the committed MitoHPC test data, covering BOTH the BAM and CRAM input paths,
# then asserts:
#   * every caller "operates": status ok AND its expected output file exists
#     (positive sample sv_del4977_h30, which carries the ~4977 bp deletion)
#   * the pipeline detects the known common deletion (>=1 caller) in the positive
#     and in the CRAM round-trip of the same sample
#   * the wild-type negative runs cleanly (specificity is reported, not gated)
###############################################################################
set -uo pipefail

IMAGE="${1:-mito-sv:ci}"
CALLERS="${2:-all}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$(mktemp -d)"
# The container runs as root, so files it writes into $OUT are root-owned; the
# in-container EXIT trap below relaxes their permissions so this host-side
# cleanup can remove them (and we tolerate any residue regardless).
trap 'rm -rf "$OUT" 2>/dev/null || sudo -n rm -rf "$OUT" 2>/dev/null || true' EXIT

POS=sv_del4977_h30          # positive: common deletion @ ~30% heteroplasmy
WT=sv_wt                    # negative: wild type, no SV
fail=0
note() { printf '\n=== %s ===\n' "$*"; }
err()  { printf 'FAIL: %s\n' "$*" >&2; fail=1; }

echo "image:   $IMAGE"
echo "callers: $CALLERS"
echo "out:     $OUT"

###############################################################################
# Run everything inside one container session (so the CRAM round-trip uses the
# image's own samtools + bundled rCRS).
###############################################################################
note "running pipeline in container"
docker run --rm \
    -v "$REPO/test/data:/data:ro" \
    -v "$OUT:/out" \
    --entrypoint bash "$IMAGE" -lc '
set -e
# Relax perms on the bind-mounted output on exit so the (non-root) host can
# clean it up — the container writes as root.
trap "chmod -R a+rwX /out 2>/dev/null || true" EXIT
RUN=/opt/pipeline/run_sample.sh
SAM="micromamba run -n mitosv samtools"

# --- BAM input path: positive + wild-type ---
"$RUN" --input /data/bams/'"$POS"'.bam --sample '"$POS"' \
       --outdir /out/'"$POS"' --threads 2 --callers '"$CALLERS"' || true
"$RUN" --input /data/bams/'"$WT"'.bam  --sample '"$WT"'  \
       --outdir /out/'"$WT"'  --threads 2 --callers '"$CALLERS"' || true

# --- CRAM input path: encode the positive against bundled rCRS, then run ---
$SAM view -b -o /out/_pos.bam /data/bams/'"$POS"'.bam
$SAM index /out/_pos.bam
$SAM view -C -T /opt/assets/rCRS.chrM.fa -o /out/'"$POS"'.cram /out/_pos.bam
$SAM index /out/'"$POS"'.cram
"$RUN" --input /out/'"$POS"'.cram --sample '"$POS"'_cram \
       --outdir /out/'"$POS"'_cram --threads 2 --callers '"$CALLERS"' || true

# --- cohort consolidation ---
micromamba run -n mitosv python /opt/pipeline/postprocess.py --root /out
' || err "container session returned non-zero"

###############################################################################
# Assertions (on the host, over the shared $OUT)
###############################################################################
note "per-caller 'operating' check (positive sample: $POS)"
declare -A EXPECT=(
    [eklipse]="eklipse/eKLIPse_deletions.csv"
    [mitosalt]="mitosalt/${POS}.mitosalt.tsv"
    [splicebreak2]="splicebreak2/${POS}_LargeMTDeletions_WGS-only_NoPositionFilter.txt"
    [mitomut]="mitomut/mitomut_results.txt"
    [mitoseek]="mitoseek/mitoseek_large_deletion.sam"
)
if [[ -f "$OUT/$POS/status.tsv" ]]; then
    echo "status.tsv:"; cat "$OUT/$POS/status.tsv"
else
    err "no status.tsv for $POS (preprocess likely failed)"
fi
for caller in eklipse mitosalt splicebreak2 mitomut mitoseek; do
    [[ "$CALLERS" == "all" || ",$CALLERS," == *",$caller,"* ]] || continue
    st="$(awk -F'\t' -v c="$caller" '$1==c{print $2}' "$OUT/$POS/status.tsv" 2>/dev/null)"
    out_ok=0; [[ -f "$OUT/$POS/${EXPECT[$caller]}" ]] && out_ok=1
    if [[ "$st" == "ok" && "$out_ok" == 1 ]]; then
        echo "  OK   $caller (status=$st, output present)"
    else
        err "$caller not operating (status='${st:-missing}', output_present=$out_ok)"
    fi
done

note "sensitivity: common (~4977 bp) deletion detection"
cd_file="$OUT/cohort_common_deletion.tsv"
[[ -f "$cd_file" ]] || err "cohort_common_deletion.tsv not produced"
if [[ -f "$cd_file" ]]; then
    pos_hits="$(awk -F'\t' -v s="$POS" '$1==s && $3==1{print $2}' "$cd_file" | sort -u)"
    cram_hits="$(awk -F'\t' -v s="${POS}_cram" '$1==s && $3==1{print $2}' "$cd_file" | sort -u)"
    echo "  $POS detected by: ${pos_hits:-<none>}"
    echo "  ${POS}_cram detected by: ${cram_hits:-<none>}"
    [[ -n "$pos_hits" ]]  || err "no caller detected the common deletion in $POS (BAM path)"
    [[ -n "$cram_hits" ]] || err "no caller detected the common deletion in ${POS}_cram (CRAM path)"
fi

note "specificity (informational): wild-type $WT"
if [[ -f "$cd_file" ]]; then
    wt_hits="$(awk -F'\t' -v s="$WT" '$1==s && $3==1{print $2}' "$cd_file" | sort -u)"
    echo "  $WT flagged common deletion by: ${wt_hits:-<none>}  (expected: <none>)"
fi

note "cohort summary"
[[ -f "$OUT/cohort_summary.txt" ]] && cat "$OUT/cohort_summary.txt"
echo; echo "caller matrix:"; [[ -f "$OUT/cohort_caller_matrix.tsv" ]] && cat "$OUT/cohort_caller_matrix.tsv"

note "result"
if [[ "$fail" == 0 ]]; then
    echo "SMOKE TEST PASSED"; exit 0
else
    echo "SMOKE TEST FAILED"; exit 1
fi
