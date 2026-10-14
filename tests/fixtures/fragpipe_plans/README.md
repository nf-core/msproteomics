# FragPipe plan fixtures

Each `<name>/` directory is the work dir (`fp_work/`) of a FragPipe 24.0 headless `--write-plan --msfragger-per-file` run (FRAGPIPEPLAN) on a saved workflow and manifest.
`<name>.workflow` is that run's saved `fragpipe.workflow`, the parser's input in `tests/python/test_plan_to_tool_configs.py`.
The data are public: PXD014415 (`LFQ-MBR-paper_perfile`), PXD017646 (`glyco-O-Pair_perfile`), PXD030670 (`glyco-N-*`), PXD047839 (`DIA_SpecLib_Quant_diaPASEF_perfile`), the CPTAC ccRCC DIA tutorial and subset (`ccRCC_DIA_tutorial`, `DIA_*`), and the `test_dda_lfq` / `test_tmtq` test data (`LFQ-MBR_test_dda_lfq`, `TMT6_test_tmtq`).

The recorded files are normalised:

- the task dir is written as `/plan`, and the run log has its licence material redacted;
- FragPipe's install is `@ROOT@/fragpipe-24.0`, with the academic build's jar names and tool versions.

`tests/python/test_plan_to_tool_configs.py` (InstallLayoutTest) checks that a real install path and jar names with an edition word give the same tool configs.
The files are byte-exact: prettier, the whitespace hooks and `nextflow lint` skip this directory.

To record a fixture, run FRAGPIPEPLAN with a FragPipe build that has `--write-plan` and `--msfragger-per-file` ([Nesvilab/FragPipe#2907](https://github.com/Nesvilab/FragPipe/pull/2907)), for example through `modules/local/fragpipeplan/tests/golden*.nf.test` with `FRAGPIPE_CONTAINER` set to that build, and copy its `fp_work/` with the normalisation above.

## Identical work dirs

`glyco-N-LFQ-MBG_perfile/` and `glyco-N-preptms_perfile/` hold the same files apart from the run log.
They come from two runs of one workflow, so the dry runs are the same.
`glyco-N-LFQ-MBG_perfile.workflow` and `glyco-N-preptms_perfile.workflow` differ in two lines only:

- the `# Pandas version` comment, which `java.util.Properties` skips;
- `workdir=`, which headless FragPipe replaces with `--workdir` (`propsFile.setProperty("workdir", Fragpipe.workdir)`, Fragpipe.java:513).

No setting differs, so no config file FragPipe writes can differ.
