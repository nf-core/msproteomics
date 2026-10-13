# EXPORTMATCHEDFRAGMENTS

Annotate every PSM in each group's `psm.tsv` with its matched fragment ions, as FragPipe's "Export matched fragments" option does.

## Description

FragPipe 24.0 runs `CmdExportMatchedFragments` when the workflow sets `tab-run.export_matched_fragments=true` (`TabRun.java:223-225,491-492`).
The command is `java -cp FP-PDV/FragAnnotator-1.0.0.jar:batmass-io-<version>.jar FragmentsExportMainClass <workDir> <threads> r`, run in the FragPipe work dir.
FragPipe orders it right after PhilosopherReport and before IonQuant, PTM-Shepherd, TMT-Integrator and SpecLibGen (`FragpipeRun.java:2411-2417,2490`).
It rewrites each `psm.tsv` in place and appends the columns `ions`, `ion_mz` and `ion_int`.
No FragPipe license is needed.
The export finds each group's `psm.tsv` by the manifest's experiment (and bioreplicate) name.
When any experiment is blank, or SpecLibGen runs, it reads only the work-dir root `psm.tsv`, so groups such as `exp_<bioreplicate>` are left unannotated, as in FragPipe.
The module reproduces this by running the jar on FragPipe's manifest.
Spectra come from the channel MSFragger searched.
A `.raw` file is listed as `.raw` when MSFragger wrote its `_uncalibrated.mzML` or `_calibrated.mzML` (native .raw reading, RC-09); the jar then reads that mzML, as in FragPipe.
Otherwise (the .raw was converted to mzML before the search, or the split-database search wrote no `_uncalibrated.mzML`) the line names the converted `<name>.mzML`.
This is a deviation: FragPipe would look for a missing `_calibrated.mzML` and the export would exit 1.

FragAnnotator exits 0 when it finds no manifest or cannot read the results ("No valid manifest file", "Error Parsing File"), for example when the manifest's blank experiment points it at a work-dir root without `psm.tsv`.
FragPipe then goes on with every `psm.tsv` unchanged, and so does the module; it prints a `WARN` line to the task's stderr.

The module builds the work-dir layout the jar reads:

| Path                                  | Content                                                                                                                       |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `fragpipe/fragpipe-files.fp-manifest` | FragPipe's work-dir manifest (`TabWorkflow.manifestSave`): `../spectra/<name>`, sanitised experiment, bioreplicate, data type |
| `fragpipe/fragpipe.workflow`          | `speclibgen.run-speclibgen=<true\|false>`, the only workflow key the export reads                                             |
| `fragpipe/<group>/psm.tsv`            | Named groups (copies of the Philosopher results)                                                                              |
| `fragpipe/psm.tsv`                    | FragPipe's blank group ""                                                                                                     |
| `fragpipe/umpire-se.params`           | Present when DIA-Umpire runs, so DIA files are read from `fragpipe/<name>_*.mzML`                                             |
| `spectra/`                            | Searched mzML files, `<name>_uncalibrated.mzML` of Bruker DDA `.d` files and `<name>_diatracer.mzML`                          |

FRAGPIPE_MATCHED_FRAGMENTS (`subworkflows/local/fragpipe_matched_fragments/`) runs the module after FRAGPIPE_INFERENCE, or after the PHILOSOPHERREPORT in FRAGPIPE_QUANT when FreeQuant or Labelquant runs.
The annotated `psm.tsv` is published in `philosopher/` in place of the report's `psm.tsv`.

## Inputs

| Channel        | Type         | Description                                                     |
| -------------- | ------------ | --------------------------------------------------------------- |
| meta           | val(map)     | Aggregate metadata map                                          |
| group_dirs     | path         | Per-group Philosopher result directories (staged as `groups/*`) |
| spectra        | path         | Spectra the manifest names (staged as `spectra/*`)              |
| umpire_spectra | path         | DIA-Umpire pseudo-spectra (staged as `fragpipe/*`)              |
| manifest       | path         | Manifest for the export                                         |
| root_group_dir | val(string)  | Stand-in directory of the blank group, `''` when there is none  |
| diaumpire      | val(boolean) | DIA-Umpire runs                                                 |

## Outputs

| Channel                | Type              | Description                                         |
| ---------------------- | ----------------- | --------------------------------------------------- |
| results_dirs           | tuple(meta, path) | `results/<group dir>/` with the annotated `psm.tsv` |
| psms                   | tuple(meta, path) | `results/<group dir>/psm.tsv`                       |
| log                    | tuple(meta, path) | `export_matched_fragments.log`                      |
| versions_fragannotator | tuple             | FragAnnotator version (topic channel)               |

## Parameters

| Setting                  | Description                                                                                                      |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------- |
| `ext.args`               | Appended after FragPipe's arguments (`<workDir> <threads> r`); empty by default, as FragPipe passes nothing else |
| `ext.fragpipe_tools_dir` | FragPipe `tools/` folder holding `FP-PDV/FragAnnotator-1.0.0.jar` and one `batmass-io-*.jar`                     |
| `--fragpipe_container`   | Image providing the FragPipe tools                                                                               |
