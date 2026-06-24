# Smoke-test summary

Functional status of each caller on the committed test data, produced by
`test/smoke_test.sh` during the CI image build. The positive control
`sv_del4977_h30` carries the ~4977 bp common deletion (m.8470_13447del).

- **ran** — the caller completed and produced its expected output file (gated: a 'no' fails the build)
- **detected common deletion** — it actually called del4977 in `sv_del4977_h30` (not gated; a miss is only a warning)

| caller | ran | detected common deletion |
|--------|:---:|:------------------------:|
| eklipse | yes | yes |
| mitosalt | yes | no |
| splicebreak2 | yes | no |
| mitomut | yes | yes |
| mitoseek | yes | yes |

Callers that ran: 5/5

## CRAM input path (`sv_del4977_h30_cram`)

Common deletion detected by: eklipse,mitomut,mitoseek

## Specificity — wild-type negative (`sv_wt`)

Common deletion flagged by: (none)  _(expected: none)_
