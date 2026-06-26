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
| splicebreak2 | yes | yes |
| mitomut | yes | yes |
| mitoseek | yes | yes |

## Caller comparison across MitoHPC scenarios

Evaluation only — how each third-party caller behaves on the diverse MitoHPC test constructs (we do not control their source, so nothing here gates the build). A **deletion-like call** = a record typed deletion/duplication with a span matching the truth event (common deletion within +/-80 bp; others within +/-250 bp). Forward-looking / ambiguous rows (origin-crossing, dup-del, sub-size) are shown but not scored.

### Accuracy (all scored scenarios)

| caller | sensitivity | specificity | precision | F1 | bal.acc | MCC | FP |
|--------|:-----------:|:-----------:|:---------:|:--:|:-------:|:---:|:--:|
| mitohpc | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |
| eklipse | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |
| mitosalt | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |
| splicebreak2 | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |
| mitomut | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |
| mitoseek | 17% (2/12) | 100% (11/11) | 100% | 0.29 | 58% | 0.30 | 0 |

### Deletions

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_del4977_h05 | del 8469–13447 · COMMON | · | · | · | · | · | · |
| sv_del4977_h30 | del 8469–13447 · COMMON | detected | detected | detected | detected | detected | detected |
| sv_del6000_h50 | del 5999–10999 | · | · | · | · | · | · |
| sv_del_13kb | del 2000–15001 | · | · | · | · | · | · |
| sv_del_45 | del 9000–9046 | · | · | · | · | · | · |
| sv_del_500 | del 8000–8501 | · | · | · | · | · | · |
| sv_dloop | del 400–6000 | · | · | · | · | · | · |
| sv_homoplasmy | del 8469–13447 · COMMON | · | · | · | · | · | · |
| sv_lowcov | del 8469–13447 · COMMON | · | · | · | · | · | · |
| sv_multidel | del 8469–13447 · COMMON | · | · | · | · | · | · |
| sv_multidel | del 5999–10999 | · | · | · | · | · | · |
| sv_del4977_h30_cram | del 8469–13447 · COMMON | detected | detected | detected | detected | detected | detected |
| spike_del4977_h20 | del 8469–13447 · COMMON | · | · | · | · | · | · |

### Controls

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_wt | wild-type (no SV) | · | · | · | · | · | · |
| NA12718 | wild-type (no SV) | · | · | · | · | · | · |
| NA12748 | wild-type (no SV) | · | · | · | · | · | · |
| NA12775 | wild-type (no SV) | · | · | · | · | · | · |

### Duplications *(evaluation-only / forward-looking)*

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_dup | dup (no deletion) | · | · | · | · | · | · |
| sv_dup_large | dup (no deletion) | · | · | · | · | · | · |

### Inversions *(evaluation-only / forward-looking)*

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_inv_large | inv (no deletion) | · | · | · | · | · | · |
| sv_inv_lowhet | inv (no deletion) | · | · | · | · | · | · |
| sv_inv_origin | inv (no deletion) | · | · | · | · | · | · |
| sv_inv_small | inv (no deletion) | · | · | · | · | · | · |

### Origin-crossing *(evaluation-only / forward-looking)*

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_del_origin_spares | delwrap 16400–100 | · | · | · | · | · | · |
| sv_origin | delwrap 16400–200 | · | · | · | · | · | · |

### Complex *(evaluation-only / forward-looking)*

| sample | truth event | mitohpc | eklipse | mitosalt | splicebreak2 | mitomut | mitoseek |
|--------|-------------|--|--|--|--|--|--|
| sv_dupdel | dupdel (no deletion) | · | · | · | · | · | · |
| sv_invdup | invdup (no deletion) | · | · | · | · | · | · |

