# PTM-Shepherd, PairScans, O-Pair and MBG cases

`cases.json` lists 80 workflows: the bundled workflows that run PTM-Shepherd or O-Pair, plus variants that set a key to a non-default value, an unformatted number or an invalid combo item, or leave it out (`capture.py case_workflow`).
`expected/<case>.json` is FragPipe 24.0's own output for each case, captured with `capture.sh` from a headless `--dry-run --export-plan` run: the shepherd.config lines and the PairScans, O-Pair and MBG argv, or FragPipe's refusal.

`plans/<case>/` is FRAGPIPEPLAN's dry-run work dir of the same case, recorded from `plans/<case>.plan_inputs.json` (`generate_specs.py`) with `RECORD_SPECS=tests/fixtures/ptmshepherd_parity/plans` and `tests/lib/record_fragpipe_plans.nf.test`, then `tests/lib/normalise_fragpipe_plan.py`.
A case FragPipe refuses has `plans/<case>.refusal.txt` instead: the task's exit status and FragPipe's error lines.
`tests/python/test_ptmshepherd_parity.py` compares the planner's ptmshepherd, pairscans, opair and mbg entries on each work dir with `expected/`, and each refusal with FragPipe's message.
