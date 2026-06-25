# Under-the-hood completion (positive control `sv_del4977_h30`)

Confirms each caller's *internal* pipeline ran to completion — not just that
the wrapper exited 0. Captured by `test/smoke_test.sh`. The full verbose logs
for every sample x caller are uploaded by CI as the `caller-diagnostics` artifact.

| caller | ran | clean exit | internal-pipeline signal |
|--------|:---:|:----------:|--------------------------|
| mitohpc | yes | yes | output present |
| eklipse | yes | yes | output present |
| mitosalt | yes | yes | native_rc=1 completed=1 split_aln=60905 paired_name_arms=60905/60905 lowscore_arms=0/60905 breakpoints=119 clusters=1 calls=0  |
| splicebreak2 | yes | no(16) | native_rc=0 junctions=1 del4977_junc=1 result_bytes=0 calls=0  |
| mitomut | yes | yes | output present |
| mitoseek | yes | no(1) | output present |

_Signal glossary — mitosalt:_ `split_aln` LAST split rows, `paired_name_arms`
arms whose query name ends /1|/2, `lowscore_arms` arms dropped by the score
filter, `breakpoints`/`clusters`/`calls` downstream survivors.
_splicebreak2:_ `junctions` MapSplice junctions, `del4977_junc` junctions
spanning ~8470..13447, `result_bytes` (0 => inner script exited before its
header), `calls` deletion rows.
