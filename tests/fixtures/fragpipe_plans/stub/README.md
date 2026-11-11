# Stub-test FragPipe plans

This directory holds the recorded FragPipe dry runs that the stub tests plan from (`--fragpipe_plan`).
In stub mode FRAGPIPEPLAN writes an empty plan, so every stub test that reaches FRAGPIPE_TOOL_CONFIGS sets `fragpipe_plan` to one of the `<name>/` directories here.

Each `<name>.plan_inputs.json` holds the inputs FRAGPIPEPLAN gets in that test: the workflow text, the manifest rows, the TMT annotation rows, the database name and the FragPipe options (`ext_args`).
`generate_specs.py` writes them from the tests, the samplesheet fixtures and the workflow files, as the pipeline builds them.
Its `sources` and `database_source` keys are informational only; the recorder does not read them.
Tests whose specs agree on everything the recorder reads share one `<name>`.

Each `<name>/` is FRAGPIPEPLAN's work dir (`fp_work`) recorded by `tests/lib/record_fragpipe_plans.nf.test` from `<name>.plan_inputs.json`, then normalised by `tests/lib/normalise_fragpipe_plan.py` (task dir written as `/plan`, licence redacted).
FRAGPIPE_TOOL_CONFIGS checks that the plan lists exactly the test's LC-MS files, with their experiment, bioreplicate, data type and FragPipe group (`checkFragpipePlanManifest`).
A spec FragPipe's dry run refuses has no `<name>/`: its stub test asserts the refusal instead.

## Re-record

Regenerate the specs after changing a test's inputs:

```bash
python3 -I tests/fixtures/fragpipe_plans/stub/generate_specs.py "$PWD"
```

Record them with a FragPipe image whose `--dry-run` writes the configs, on AWS Batch (not locally):

```bash
FRAGPIPE_CONTAINER=<image> nf-test test tests/lib/record_fragpipe_plans.nf.test
```

Then normalise each task's `fp_work` (the task's name, `FRAGPIPEPLAN (<name>)`, is in its `.command.run`) into its fixture directory:

```bash
python3 tests/lib/normalise_fragpipe_plan.py <task dir>/fp_work tests/fixtures/fragpipe_plans/stub/<name>
```

The normaliser stops when a file still holds the licence.
Before committing, check that no fixture file holds an internal path (buckets, image names, host or scratch directories).
