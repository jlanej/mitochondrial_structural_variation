#!/usr/bin/env bash
###############################################################################
# smoke_test.sh — robust functional CI for the built image.
#
#   bash test/smoke_test.sh [IMAGE] [CALLERS]
#
# Runs the full per-sample pipeline (preprocess + all callers) inside the image
# on the committed MitoHPC test data, covering BOTH the BAM and CRAM input paths.
#
# HARD FAILS (block the build) — "does it run?" is the bar:
#   * preprocess failed (no status.tsv), OR
#   * post-processing produced no cohort tables, OR
#   * ANY caller did not run (status != ok or no output file); its log is printed.
# LOUD WARNINGS (do NOT block the build) — "did it call?":
#   * the known common deletion was not detected by a caller that ran.
# Also writes deterministic example outputs + a SMOKE_SUMMARY.md (operated vs
# detected per caller) to test/example_output/ for CI to commit.
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
#
# Gating policy:
#   HARD FAIL  — a caller did not RUN (status != ok or no output file), or
#                preprocess/post-processing broke. "It runs" is the bar.
#   WARN only  — a caller ran but did not DETECT the known common deletion.
# Logs of any non-running caller are printed so they can be diagnosed.
###############################################################################
warns=0
warn() { printf 'WARNING: %s\n' "$*" >&2; warns=$((warns+1)); }

note "per-caller 'did it run' check (positive sample: $POS)"
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
operating=0
declare -A OP=()
for caller in eklipse mitosalt splicebreak2 mitomut mitoseek; do
    if [[ "$CALLERS" != "all" && ",$CALLERS," != *",$caller,"* ]]; then OP[$caller]=skip; continue; fi
    st="$(awk -F'\t' -v c="$caller" '$1==c{print $2}' "$OUT/$POS/status.tsv" 2>/dev/null)"
    out_ok=0; [[ -f "$OUT/$POS/${EXPECT[$caller]}" ]] && out_ok=1
    if [[ "$st" == "ok" && "$out_ok" == 1 ]]; then
        echo "  OK   $caller ran (status=$st, output present)"; operating=$((operating+1)); OP[$caller]=yes
    else
        OP[$caller]=no
        err "$caller did NOT run (status='${st:-missing}', output_present=$out_ok)"
        clog="$OUT/$POS/$caller/$caller.log"
        if [[ -f "$clog" ]]; then
            echo "----- $caller log (tail) -----------------------------------------" >&2
            tail -n 30 "$clog" >&2
            echo "------------------------------------------------------------------" >&2
        fi
    fi
done

note "sensitivity: common (~4977 bp) deletion detection (WARN only, not gated)"
cd_file="$OUT/cohort_common_deletion.tsv"
[[ -f "$cd_file" ]] || err "cohort_common_deletion.tsv not produced (post-processing broken)"
if [[ -f "$cd_file" ]]; then
    pos_hits="$(awk -F'\t' -v s="$POS" '$1==s && $3==1{print $2}' "$cd_file" | sort -u | paste -sd, -)"
    cram_hits="$(awk -F'\t' -v s="${POS}_cram" '$1==s && $3==1{print $2}' "$cd_file" | sort -u | paste -sd, -)"
    echo "  $POS detected by: ${pos_hits:-<none>}"
    echo "  ${POS}_cram detected by: ${cram_hits:-<none>}"
    [[ -n "$pos_hits" ]]  || warn "no caller detected the common deletion in $POS (BAM path)"
    [[ -n "$cram_hits" ]] || warn "no caller detected the common deletion in ${POS}_cram (CRAM path)"
fi

note "specificity (informational): wild-type $WT"
if [[ -f "$cd_file" ]]; then
    wt_hits="$(awk -F'\t' -v s="$WT" '$1==s && $3==1{print $2}' "$cd_file" | sort -u | paste -sd, -)"
    echo "  $WT flagged common deletion by: ${wt_hits:-<none>}  (expected: <none>)"
fi

note "cohort summary"
[[ -f "$OUT/cohort_summary.txt" ]] && cat "$OUT/cohort_summary.txt"
echo; echo "caller matrix:"; [[ -f "$OUT/cohort_caller_matrix.tsv" ]] && cat "$OUT/cohort_caller_matrix.tsv"

###############################################################################
# Save deterministic example outputs into the repo (CI commits these on a
# successful build). We copy only small text result files — NOT logs/status
# (timestamps churn) and NOT BAM/FASTQ intermediates.
###############################################################################
note "saving example outputs -> test/example_output"
EXDIR="$REPO/test/example_output"
rm -rf "$EXDIR"; mkdir -p "$EXDIR"
cp "$OUT"/cohort_sv_calls.tsv "$OUT"/cohort_common_deletion.tsv \
   "$OUT"/cohort_caller_matrix.tsv "$OUT"/cohort_summary.txt "$EXDIR"/ 2>/dev/null || true
for s in "$POS" "$WT" "${POS}_cram"; do
    sd="$OUT/$s"; [[ -d "$sd" ]] || continue
    for caller in eklipse mitosalt splicebreak2 mitomut mitoseek; do
        cdir="$sd/$caller"; [[ -d "$cdir" ]] || continue
        dest="$EXDIR/$s/$caller"; mkdir -p "$dest"
        # only top-level result text files (skip work/, sb_install/, *.bam,
        # input lists, etc.)
        find "$cdir" -maxdepth 1 -type f ! -name 'bam_list.tsv' \
            \( -name '*.csv' -o -name '*.tsv' -o -name '*.txt' -o -name '*.sam' \) \
            -exec cp {} "$dest/" \; 2>/dev/null || true
        rmdir "$dest" 2>/dev/null || true   # drop empty dirs (failed callers)
    done
    rmdir "$EXDIR/$s" 2>/dev/null || true
done

# Human-readable "what is working vs just running" summary (deterministic — no
# timestamps — so it only changes when caller behaviour changes).
detected_by() {  # sample caller -> "yes"/"no"
    awk -F'\t' -v s="$1" -v c="$2" '$1==s && $2==c && $3==1{f=1} END{print (f?"yes":"no")}' \
        "$cd_file" 2>/dev/null || echo "no"
}
summary_md="$EXDIR/SMOKE_SUMMARY.md"
{
    echo "# Smoke-test summary"
    echo
    echo "Functional status of each caller on the committed test data, produced by"
    echo '`test/smoke_test.sh` during the CI image build. The positive control'
    echo "\`$POS\` carries the ~4977 bp common deletion (m.8470_13447del)."
    echo
    echo "- **ran** — the caller completed and produced its expected output file (gated: a 'no' fails the build)"
    echo "- **detected common deletion** — it actually called del4977 in \`$POS\` (not gated; a miss is only a warning)"
    echo
    echo "| caller | ran | detected common deletion |"
    echo "|--------|:---:|:------------------------:|"
    for caller in eklipse mitosalt splicebreak2 mitomut mitoseek; do
        echo "| $caller | ${OP[$caller]:-?} | $(detected_by "$POS" "$caller") |"
    done
    echo
    echo "Callers that ran: ${operating}/5"
    echo
    echo "## CRAM input path (\`${POS}_cram\`)"
    echo
    echo "Common deletion detected by: ${cram_hits:-(none)}"
    echo
    echo "## Specificity — wild-type negative (\`$WT\`)"
    echo
    echo "Common deletion flagged by: ${wt_hits:-(none)}  _(expected: none)_"
} > "$summary_md"
echo "--- SMOKE_SUMMARY.md ---"; cat "$summary_md"
echo "--- example_output files ---"
( cd "$EXDIR" && find . -type f | sort | sed 's#^\./#  #' ) || true

note "result"
echo "operating callers on $POS: ${operating:-0}/5    warnings: $warns"
if [[ "$fail" == 0 ]]; then
    echo "SMOKE TEST PASSED"; exit 0
else
    echo "SMOKE TEST FAILED"; exit 1
fi
