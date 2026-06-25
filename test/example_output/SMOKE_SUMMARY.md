# Smoke-test summary

Functional status of each caller across the MitoHPC test cohort, produced by
`test/smoke_test.sh` during the CI image build (scope: full).

## Operating + common-deletion detection (positive control `sv_del4977_h30`)

- **ran** — completed and produced its expected output file (GATED — must pass)
- **detected common deletion** — called del4977 in `sv_del4977_h30` (evaluation only)

| caller | ran | detected common deletion |
|--------|:---:|:------------------------:|
| mitohpc | yes | yes |
| eklipse | yes | yes |
| mitosalt | yes | yes |
| splicebreak2 | yes | no |
| mitomut | yes | yes |
| mitoseek | yes | yes |

## Caller comparison across MitoHPC scenarios

Evaluation only — a record of how each third-party caller behaves on the diverse test constructs (we do not control their source, so nothing here gates the build). **detected** = at least one caller matched the truth deletion (common deletion within +/-80 bp; others within +/-250 bp).

| sample | truth event | detected | callers |
|--------|-------------|:--------:|---------|
| sv_del4977_h05 | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_del4977_h30 | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_del6000_h50 | del 5999-10999 | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_dloop | del 400-6000 | yes | eklipse, mitohpc, mitomut, mitoseek |
| sv_dup | dup 6000-7000 | yes | mitohpc, mitomut, mitosalt, mitoseek |
| sv_homoplasmy | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_lowcov | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt |
| sv_multidel | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_multidel | del 5999-10999 | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| sv_origin | del 16400-200 | no | (none) |
| sv_wt | wild-type (no SV) | - | n/a (specificity sample) |
| sv_del4977_h30_cram | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| spike_del4977_h20 | del 8469-13447 [COMMON] | yes | eklipse, mitohpc, mitomut, mitosalt, mitoseek |
| NA12718 | wild-type (no SV) | - | n/a (specificity sample) |

**Sensitivity — truth events no caller detected (observation):**

- sv_origin del 16400-200 detected by no caller

