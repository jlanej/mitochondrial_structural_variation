#!/usr/bin/env bash
###############################################################################
# smoke_test.sh — functional CI mirroring the MitoHPC sv-calling scenarios.
#
#   bash test/smoke_test.sh [IMAGE] [CALLERS] [SCOPE]
#     IMAGE   : container image (default mito-sv:ci)
#     CALLERS : caller list (default all)
#     SCOPE   : full (default) = all 10 mock BAMs + CRAM round-trip + real BAMs
#                                + degenerate-input robustness;
#               quick          = sv_del4977_h30 + sv_wt + CRAM only
#
# Runs the full per-sample pipeline inside the image over the committed MitoHPC
# test cohort (the same diverse constructs MitoHPC's own caller is tested on:
# common deletion at varying VAF/depth, non-repeat deletion, D-loop deletion,
# multi-deletion, tandem duplication, origin-crossing deletion, wild-type, plus
# real 1000G + a del4977 spike-in), then asserts:
#
#   HARD FAILS (block the build) — only things the PIPELINE controls:
#     * any caller did not RUN (produce output) on the canonical positive
#       sv_del4977_h30  ("are the callers operating?")
#     * post-processing produced no cohort tables
#     * a degenerate input (wrong-contig / empty BAM) did not fail cleanly
#   EVALUATION ONLY (never blocks): how each third-party caller DETECTS the
#     diverse constructs is recorded as a caller-comparison matrix. We do not
#     control the callers' source, so their sensitivity/specificity behaviour is
#     reported (in SMOKE_SUMMARY.md), not asserted.
#
# Writes test/example_output/ (deterministic result files + SMOKE_SUMMARY.md +
# any FAILED caller logs) for CI to commit.
###############################################################################
set -uo pipefail

IMAGE="${1:-mito-sv:ci}"
CALLERS="${2:-all}"
SCOPE="${3:-${MITO_SV_SMOKE_SCOPE:-full}}"
# SUITE = which mock BAMs the scenario cohort runs. "all" = the full MitoHPC SV
# suite (deletions + duplications + inversions + complex + origin); "del" =
# deletions + controls only (skips the forward-looking dup/inv/complex BAMs) — an
# easy lever to shorten CI if the full suite runs too long. Flip via the 4th arg
# or MITO_SV_SUITE.
SUITE="${4:-${MITO_SV_SUITE:-all}}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

POS=sv_del4977_h30          # canonical positive control (common deletion @30%)
WT=sv_wt                    # wild-type negative
TRUTH="$REPO/test/data/truth.tsv"
REALS=(spike_del4977_h20 NA12718 NA12748 NA12775)   # committed real 1000G + spike
fail=0; warns=0
note() { printf '\n=== %s ===\n' "$*"; }
err()  { printf 'FAIL: %s\n' "$*" >&2; fail=1; }
warn() { printf 'WARNING: %s\n' "$*" >&2; warns=$((warns+1)); }

###############################################################################
# Phase: full (default, one job) | shard (run a sample subset) | consolidate
# (merge shard outputs). Matrix CI fans the cohort out across runner jobs:
#   shard       MITO_SV_RUN_ONLY="s1 s2 …" MITO_SV_OUT=DIR [MITO_SV_RUN_EXTRAS=1]
#               -> run ONLY those samples into DIR, skip post-processing/gates.
#   consolidate MITO_SV_MERGED_OUT=DIR
#               -> skip the container; post-process + gate + report the merged DIR.
# Per-sample runs are independent (separate outdirs); sharding across JOBS (not
# in-process on one runner) parallelizes safely without OOMing the heavy callers.
###############################################################################
SHARD="${MITO_SV_RUN_ONLY:-}"
MERGED="${MITO_SV_MERGED_OUT:-}"
EXTRAS="${MITO_SV_RUN_EXTRAS:-}"
if [[ -n "$MERGED" ]]; then
    PHASE=consolidate; OUT="$MERGED"
elif [[ -n "$SHARD" ]]; then
    PHASE=shard; OUT="${MITO_SV_OUT:?set MITO_SV_OUT for shard mode}"
else
    PHASE=full; OUT="$(mktemp -d)"
    trap 'rm -rf "$OUT" 2>/dev/null || sudo -n rm -rf "$OUT" 2>/dev/null || true' EXIT
fi
mkdir -p "$OUT"

# Mock BAM samples for the selected suite (del = only samples carrying a
# deletion / wild-type event; all = every committed mock BAM).
select_bams() {
    local b s
    for b in "$REPO"/test/data/bams/*.bam; do
        s="$(basename "$b" .bam)"
        if [[ "$SUITE" == del ]]; then
            awk -v s="$s" '$1==s && ($2=="del"||$2=="delwrap"||$2=="none"){f=1}
                           END{exit !f}' "$TRUTH" && echo "$s"
        else
            echo "$s"
        fi
    done
}
MOCK_SAMPLES=(); while IFS= read -r _s; do MOCK_SAMPLES+=("$_s"); done < <(select_bams)

# The full sample set this run evaluates (host-side), and which BAMs to RUN now.
if [[ "$SCOPE" == quick ]]; then
    SAMPLES=("$POS" "$WT" "${POS}_cram")
else
    SAMPLES=("${MOCK_SAMPLES[@]}" "${POS}_cram" "${REALS[@]}")
fi
if [[ "$PHASE" == shard ]]; then
    read -ra RUN_SAMPLES <<< "$SHARD"
    DO_CRAM=0; DO_DEGEN=0
    [[ -n "$EXTRAS" ]] && { DO_CRAM=1; DO_DEGEN=1; }
elif [[ "$SCOPE" == quick ]]; then
    RUN_SAMPLES=("$POS" "$WT"); DO_CRAM=1; DO_DEGEN=0
else
    RUN_SAMPLES=("${MOCK_SAMPLES[@]}" "${REALS[@]}"); DO_CRAM=1; DO_DEGEN=1
fi

echo "image:$IMAGE  callers:$CALLERS  scope:$SCOPE  suite:$SUITE  phase:$PHASE"
[[ "$PHASE" != consolidate ]] && echo "run (${#RUN_SAMPLES[@]} BAMs, cram=$DO_CRAM degen=$DO_DEGEN): ${RUN_SAMPLES[*]}"
echo "out: $OUT"

###############################################################################
# Run the requested samples inside the image (skipped in consolidate phase).
###############################################################################
if [[ "$PHASE" != consolidate ]]; then
note "running samples in container (phase=$PHASE)"
docker run --rm \
    -v "$REPO/test/data:/data:ro" \
    -v "$OUT:/out" \
    -e CALLERS="$CALLERS" -e POS="$POS" -e DO_CRAM="$DO_CRAM" -e DO_DEGEN="$DO_DEGEN" \
    -e RUN_SAMPLES="${RUN_SAMPLES[*]}" \
    --entrypoint bash "$IMAGE" -lc '
set -e
trap "chmod -R a+rwX /out 2>/dev/null || true" EXIT
RUN=/opt/pipeline/run_sample.sh
SAM="micromamba run -n mitosv samtools"
run_one() { "$RUN" --input "$1" --sample "$2" --outdir /out/"$2" --threads 2 --callers "$CALLERS" || true; }
resolve() {  # sample name -> input BAM path (mock vs committed real)
    case "$1" in
        spike_del4977_h20|NA12718|NA12748|NA12775) echo /data/real/"$1".chrM.bam;;
        *) echo /data/bams/"$1".bam;;
    esac
}
for s in $RUN_SAMPLES; do run_one "$(resolve "$s")" "$s"; done

if [ "$DO_CRAM" = 1 ]; then
    # CRAM round-trip of the canonical positive (decoded offline via the seeded
    # rCRS reference cache).
    $SAM view -b -o /out/_pos.bam /data/bams/"$POS".bam; $SAM index /out/_pos.bam
    $SAM view -C -T /opt/assets/rCRS.chrM.fa -o /out/"$POS".cram /out/_pos.bam; $SAM index /out/"$POS".cram
    run_one /out/"$POS".cram "${POS}_cram"
fi

if [ "$DO_DEGEN" = 1 ]; then
    # Degenerate inputs (robustness): must fail cleanly, never hang/traceback.
    : > /out/_degen.txt
    $SAM view -h /data/bams/sv_wt.bam | sed "s/SN:chrM/SN:chrZ/; s/\tchrM\t/\tchrZ\t/" \
        | $SAM view -b -o /out/_wrongcontig.bam -
    if "$RUN" --input /out/_wrongcontig.bam --sample degen_wrongcontig \
              --outdir /out/_degen_wc --threads 2 --callers eklipse > /out/_degen_wc.log 2>&1; then rc=0; else rc=$?; fi
    echo "wrongcontig $rc" >> /out/_degen.txt
    $SAM view -H /data/bams/sv_wt.bam | $SAM view -b -o /out/_empty.bam -
    if "$RUN" --input /out/_empty.bam --sample degen_empty \
              --outdir /out/_degen_empty --threads 2 --callers eklipse > /out/_degen_empty.log 2>&1; then rc=0; else rc=$?; fi
    echo "empty $rc" >> /out/_degen.txt
fi
' || err "container session returned non-zero"
fi

# Shard phase: outputs are in $OUT for artifact upload — stop before consolidation.
if [[ "$PHASE" == shard ]]; then
    note "shard complete"
    echo "$(find "$OUT" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ') sample dir(s) in $OUT"
    [[ "$fail" == 0 ]] && exit 0 || exit 1
fi

###############################################################################
# Post-process the (possibly merged) cohort — host-side, pure stdlib.
###############################################################################
note "post-processing cohort -> $OUT"
python3 "$REPO/pipeline/postprocess.py" --root "$OUT" || err "postprocess failed"

###############################################################################
# 1. Operating gate (HARD): every caller ran on the canonical positive.
###############################################################################
note "per-caller 'did it run' check (positive sample: $POS)"
declare -A EXPECT=(
    [mitohpc]="mitohpc/mitohpc.sv.tab"
    [eklipse]="eklipse/eKLIPse_deletions.csv"
    [mitosalt]="mitosalt/${POS}.mitosalt.tsv"
    [splicebreak2]="splicebreak2/${POS}_LargeMTDeletions_WGS-only_NoPositionFilter.txt"
    [mitomut]="mitomut/mitomut_results.txt"
    [mitoseek]="mitoseek/mitoseek_large_deletion.sam"
)
[[ -f "$OUT/$POS/status.tsv" ]] && { echo "status.tsv:"; cat "$OUT/$POS/status.tsv"; } \
    || err "no status.tsv for $POS (preprocess likely failed)"
# Format a caller's runstatus sidecar (caller<tab>k=v ...) into "k=v k=v".
runstat_kv() {
    [[ -f "$1" ]] || { echo "(no sidecar)"; return; }
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[i]=$i} NR==2{s="";for(i=2;i<=NF;i++)s=s h[i]"="$i" "; print s}' "$1"
}
declare -A OP=()
for caller in mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek; do
    if [[ "$CALLERS" != "all" && ",$CALLERS," != *",$caller,"* ]]; then OP[$caller]=skip; continue; fi
    st="$(awk -F'\t' -v c="$caller" '$1==c{print $2}' "$OUT/$POS/status.tsv" 2>/dev/null)"
    out_ok=0; [[ -f "$OUT/$POS/${EXPECT[$caller]}" ]] && out_ok=1
    rs="$OUT/$POS/$caller/$caller.runstatus"
    if [[ "$st" == "ok" && "$out_ok" == 1 ]]; then
        echo "  OK   $caller ran"; OP[$caller]=yes
    else
        OP[$caller]=no
        err "$caller did NOT run (status='${st:-missing}', output_present=$out_ok)"
        clog="$OUT/$POS/$caller/$caller.log"
        [[ -f "$clog" ]] && { echo "----- $caller log (tail) -----" >&2; tail -n 30 "$clog" >&2; }
    fi
    # Surface the internal-pipeline signal so "ran" vs "ran but did no real work"
    # is visible at a glance (the diagnostics artifact has the full logs).
    [[ -f "$rs" ]] && echo "       under-the-hood: $(runstat_kv "$rs")"
done

###############################################################################
# 2. Caller comparison across scenarios (EVALUATION ONLY — never gates). Builds
#    the scenario x caller matrix from the cohort table on the host.
###############################################################################
note "caller comparison across MitoHPC scenarios (evaluation, not gated)"
cohort="$OUT/cohort_sv_calls.tsv"
scen_md="$OUT/scenario_matrix.md"
[[ -f "$OUT/cohort_common_deletion.tsv" && -f "$cohort" ]] \
    || err "post-processing produced no cohort tables"
if [[ -f "$cohort" ]]; then
    python3 "$REPO/test/check_scenarios.py" --calls "$cohort" \
        --truth "$REPO/test/data/truth.tsv" --out-md "$scen_md" \
        --samples "${SAMPLES[*]}" || true
fi

###############################################################################
# 3. Degenerate-input robustness (HARD): wrong-contig / empty BAM must fail
#    cleanly (non-zero exit, no Python traceback).
###############################################################################
if [[ "$SCOPE" != quick ]]; then
    note "degenerate-input robustness"
    if [[ -f "$OUT/_degen.txt" ]]; then
        cat "$OUT/_degen.txt"
        while read -r name rc; do
            log="$OUT/_degen_${name/wrongcontig/wc}.log"
            [[ "$name" == empty ]] && log="$OUT/_degen_empty.log"
            if [[ "$rc" == 0 ]]; then
                err "degenerate input '$name' did NOT fail (exit 0)"
            elif grep -q "Traceback (most recent call last)" "$log" 2>/dev/null; then
                err "degenerate input '$name' failed with a Python traceback (not clean)"
                tail -n 15 "$log" >&2
            else
                echo "  OK   '$name' failed cleanly (exit $rc, no traceback)"
            fi
        done < "$OUT/_degen.txt"
    else
        err "degenerate-input results missing"
    fi
fi

###############################################################################
# Cohort summary to console.
###############################################################################
note "cohort summary"; [[ -f "$OUT/cohort_summary.txt" ]] && cat "$OUT/cohort_summary.txt"

###############################################################################
# Save example outputs + SMOKE_SUMMARY.md (CI commits these).
###############################################################################
note "saving example outputs -> test/example_output"
EXDIR="$REPO/test/example_output"; rm -rf "$EXDIR"; mkdir -p "$EXDIR"
cp "$OUT"/cohort_sv_calls.tsv "$OUT"/cohort_common_deletion.tsv \
   "$OUT"/cohort_caller_matrix.tsv "$OUT"/cohort_summary.txt "$EXDIR"/ 2>/dev/null || true
for s in "${SAMPLES[@]}"; do
    sd="$OUT/$s"; [[ -d "$sd" ]] || continue
    for caller in mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek; do
        cdir="$sd/$caller"; [[ -d "$cdir" ]] || continue
        dest="$EXDIR/$s/$caller"; mkdir -p "$dest"
        find "$cdir" -maxdepth 1 -type f ! -name 'bam_list.tsv' \
            \( -name '*.csv' -o -name '*.tsv' -o -name '*.txt' -o -name '*.sam' \) \
            -exec cp {} "$dest/" \; 2>/dev/null || true
        if [ -z "$(ls -A "$dest" 2>/dev/null)" ] && [ -f "$cdir/$caller.log" ]; then
            cp "$cdir/$caller.log" "$dest/${caller}.FAILED.log"
        fi
        rmdir "$dest" 2>/dev/null || true
    done
    rmdir "$EXDIR/$s" 2>/dev/null || true
done

# SMOKE_SUMMARY.md: operating table + scenario matrix.
summary_md="$EXDIR/SMOKE_SUMMARY.md"
detected_by() {  # sample caller -> yes/no (common deletion)
    awk -F'\t' -v s="$1" -v c="$2" '$1==s && $2==c && $3==1{f=1} END{print (f?"yes":"no")}' \
        "$OUT/cohort_common_deletion.tsv" 2>/dev/null || echo "no"
}
{
    echo "# Smoke-test summary"
    echo
    echo "Functional status of each caller across the MitoHPC test cohort, produced by"
    echo '`test/smoke_test.sh` during the CI image build (scope: '"$SCOPE"').'
    echo
    echo "## Operating + common-deletion detection (positive control \`$POS\`)"
    echo
    echo "- **ran** — completed and produced its expected output file (GATED — must pass)"
    echo "- **detected common deletion** — called del4977 in \`$POS\` (evaluation only)"
    echo
    echo "| caller | ran | detected common deletion |"
    echo "|--------|:---:|:------------------------:|"
    for caller in mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek; do
        echo "| $caller | ${OP[$caller]:-?} | $(detected_by "$POS" "$caller") |"
    done
    echo
    [[ -f "$scen_md" ]] && cat "$scen_md"
} > "$summary_md"
echo "--- SMOKE_SUMMARY.md ---"; cat "$summary_md"

###############################################################################
# Under-the-hood completion record (committed, low-churn) + full verbose
# diagnostics bundle (CI artifact). Lets us confirm each caller's INTERNAL
# pipeline ran — not just that the wrapper exited 0 — and see where a caller
# that "ran" produced no call.
###############################################################################
note "under-the-hood completion (positive sample) + diagnostics bundle"
uth_md="$EXDIR/UNDER_THE_HOOD.md"
{
    echo "# Under-the-hood completion (positive control \`$POS\`)"
    echo
    echo "Confirms each caller's *internal* pipeline ran to completion — not just that"
    echo 'the wrapper exited 0. Captured by `test/smoke_test.sh`. The full verbose logs'
    echo 'for every sample x caller are uploaded by CI as the `caller-diagnostics` artifact.'
    echo
    echo "| caller | ran | clean exit | internal-pipeline signal |"
    echo "|--------|:---:|:----------:|--------------------------|"
    for caller in mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek; do
        rs="$OUT/$POS/$caller/$caller.runstatus"
        clog="$OUT/$POS/$caller/$caller.log"
        # "clean exit" = wrapper log has no fatal error markers (traceback /
        # command-not-found / unhandled exception). Verifies the tool didn't
        # error its way to a 0-exit.
        nerr=0
        if [[ -f "$clog" ]]; then
            nerr="$(grep -ciE 'traceback \(most recent call last\)|command not found|no such file or directory|exception in thread|core dumped' "$clog" 2>/dev/null)" || nerr=0
        fi
        clean=$([[ "${nerr:-0}" -eq 0 ]] && echo yes || echo "no($nerr)")
        if [[ -f "$rs" ]]; then sig="$(runstat_kv "$rs")"
        else sig="output $([[ -f "$OUT/$POS/${EXPECT[$caller]:-_none_}" ]] && echo present || echo MISSING)"; fi
        echo "| $caller | ${OP[$caller]:-?} | $clean | ${sig} |"
    done
    echo
    echo "_Signal glossary — mitosalt:_ \`split_aln\` LAST split rows, \`paired_name_arms\`"
    echo "arms whose query name ends /1|/2, \`lowscore_arms\` arms dropped by the score"
    echo "filter, \`breakpoints\`/\`clusters\`/\`calls\` downstream survivors."
    echo "_splicebreak2:_ \`junctions\` MapSplice junctions, \`del4977_junc\` junctions"
    echo "spanning ~8470..13447, \`result_bytes\` (0 => inner script exited before its"
    echo "header), \`calls\` deletion rows."
} > "$uth_md"
echo "--- UNDER_THE_HOOD.md ---"; cat "$uth_md"

# TEMPORARY debug: commit the positive sample's full Splice-Break2 wrapper log
# (its under-the-hood stage map + the inner Splice-Break2_0725.sh nohup tails) so
# the inner-script failure point is visible IN-REPO (the verbose bundle below is
# artifact-only). Remove once Splice-Break2 produces non-empty output.
sb_log="$OUT/$POS/splicebreak2/splicebreak2.log"
[[ -f "$sb_log" ]] && cp "$sb_log" "$EXDIR/SPLICEBREAK2_DEBUG.log" && echo "wrote SPLICEBREAK2_DEBUG.log"

# Verbose bundle: every caller's full log + sidecar + small native logs/
# intermediates (no BAM/FASTQ/bigwig). Gitignored; uploaded as a CI artifact.
DIAG="$REPO/test/_diagnostics"; rm -rf "$DIAG"; mkdir -p "$DIAG"
for s in "${SAMPLES[@]}"; do
    sd="$OUT/$s"; [[ -d "$sd" ]] || continue
    [[ -f "$sd/status.tsv" ]] && { mkdir -p "$DIAG/$s"; cp "$sd/status.tsv" "$DIAG/$s/" 2>/dev/null || true; }
    for caller in mitohpc eklipse mitosalt splicebreak2 mitomut mitoseek; do
        cdir="$sd/$caller"; [[ -d "$cdir" ]] || continue
        dest="$DIAG/$s/$caller"; mkdir -p "$dest"
        while IFS= read -r f; do
            rel="${f#"$cdir"/}"; mkdir -p "$dest/$(dirname "$rel")"; cp "$f" "$dest/$rel" 2>/dev/null || true
        done < <(find "$cdir" -type f \( -name '*.log' -o -name '*.runstatus' -o -name '*.Rout' \
                    -o -name '*.breakpoint' -o -name '*.cluster' -o -name '*junction*.txt' \
                    -o -name '*_nohup.log' -o -name '*LargeMTDeletions*.txt' \) -size -3M 2>/dev/null)
        [[ -z "$(ls -A "$dest" 2>/dev/null)" ]] && rmdir "$dest" 2>/dev/null || true
    done
    [[ -z "$(ls -A "$DIAG/$s" 2>/dev/null)" ]] && rmdir "$DIAG/$s" 2>/dev/null || true
done
( cd "$DIAG" && find . -type f | sort > _MANIFEST.txt ) 2>/dev/null || true
echo "diagnostics bundle: $(find "$DIAG" -type f 2>/dev/null | wc -l | tr -d ' ') files -> test/_diagnostics"

###############################################################################
# Interactive comparison report -> docs/index.html (CI commits it).
###############################################################################
note "generating docs/index.html (interactive caller comparison)"
if [[ -f "$cohort" && -f "$OUT/cohort_runtime.tsv" ]]; then
    python3 "$REPO/pipeline/make_report.py" \
        --calls "$cohort" --runtime "$OUT/cohort_runtime.tsv" \
        --truth "$REPO/test/data/truth.tsv" --scenarios "$REPO/test/data/scenarios.json" \
        --samples "${SAMPLES[*]}" \
        --out "$REPO/docs/index.html" --scope "$SCOPE" \
        --image "$IMAGE" --generated "${MITO_SV_GENERATED:-$(date -u +%Y-%m-%d)}" \
        && echo "wrote docs/index.html" || warn "report generation failed"
else
    warn "cohort tables missing — skipping report"
fi

note "result"
echo "warnings: $warns"
if [[ "$fail" == 0 ]]; then echo "SMOKE TEST PASSED"; exit 0; else echo "SMOKE TEST FAILED"; exit 1; fi
