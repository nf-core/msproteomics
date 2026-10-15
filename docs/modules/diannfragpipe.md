# DIANNFRAGPIPE

FragPipe's DIA-Quant step: the DIA-NN binary bundled with FragPipe 24.0 with the arguments of `CmdDiann`, and FragPipe's post-steps.

## Description

The binary is `tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8` (DIA-NN 1.8.2 Beta 8); the module stops when it reports another version.
`LD_PRELOAD` is resolved as FragPipe does (`diann_so/libm.so.6`, `libstdc++.so.6`, then `libgomp.so.1.0.0` when the bare binary fails).
The inputs go to `filelist_diann.txt` as `--f <path>` lines, largest first, and DIA-NN runs with `--cfg filelist_diann.txt--`.
After DIA-NN the `.quant` files are deleted, the `.speclib` is moved to `dia-quant-output/report-tsv.speclib`, and FragPipe's `Propagation`, `DiannToMsstats` (MSstats on) and `FragReporter` run; plexDIA runs `PlexDiaHelper` first.

## Inputs

| Channel               | Type                             | Description                                                  |
| --------------------- | -------------------------------- | ------------------------------------------------------------ |
| meta                  | val(map)                         | Group metadata                                               |
| library               | path (stageAs `library_input/*`) | SpecLibGen's `library.tsv` or the workflow's `diann.library` |
| spectra               | path (stageAs `spectra/*`)       | DIA and DIA-Quant LC-MS files                                |
| psm                   | path (stageAs `psm.tsv`)         | The group's psm.tsv                                          |
| manifest              | path                             | `fragpipe-files.fp-manifest`                                 |
| experiment_annotation | val(string)                      | `experiment_annotation.tsv` content                          |
| diann_config          | val(string)                      | JSON of the parsed `diann` entry                             |
| fasta                 | path                             | Protein database (`--relaxed-prot-inf --fasta`)              |

## Outputs

| Channel               | Type             | Description                                         |
| --------------------- | ---------------- | --------------------------------------------------- |
| results_files         | tuple(val, path) | every file of `dia-quant-output/`                   |
| report                | tuple(val, path) | `report.tsv`                                        |
| pr_matrix             | tuple(val, path) | `report.pr_matrix.tsv`                              |
| pg_matrix             | tuple(val, path) | `report.pg_matrix.tsv`                              |
| speclib               | tuple(val, path) | `report-tsv.speclib` (optional)                     |
| msstats               | tuple(val, path) | `msstats.csv` (optional)                            |
| experiment_annotation | tuple(val, path) | `experiment_annotation.tsv`                         |
| filelist              | tuple(val, path) | `filelist_diann.txt`                                |
| argv                  | tuple(val, path) | `diann_argv.txt`, the command one argument per line |
| versions_diann        | topic: versions  | DIA-NN version                                      |
| versions_fragpipe     | topic: versions  | FragPipe version                                    |

## Parameters

| Parameter              | Default                                           | Description                           |
| ---------------------- | ------------------------------------------------- | ------------------------------------- |
| ext.args               | ''                                                | Extra DIA-NN arguments before `--cfg` |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory              |

## Container

The FragPipe image given with `--fragpipe_container`.

## References

- [DIA-NN](https://github.com/vdemichev/DiaNN)
- [FragPipe homepage](https://fragpipe.nesvilab.org/)

## Notes

- `--no-quant-files` is not passed: the FragPipe 24.0 release does not pass it (its run log has no such flag), although the `24.0` source tag adds it after `--smart-profiling`; it was removed from the source in commit `91bb5c92f`.
- Inputs of equal size keep the LC-MS table order in `filelist_diann.txt`; FragPipe keeps them in `Set` iteration order, which is not reproducible.
- A workflow with transfer learning (`transfer-learning.run-transfer-learning=true`) stops before the run (`PARSEFRAGPIPEWORKFLOW`): the pipeline has no TransferLearning step, which needs an external prediction server, so DIA-NN never reads its predicted library `fragpipe-predicted-speclib.<format>` (CmdDiann.java:131-135, :321-322).
