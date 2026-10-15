# SPECLIBGEN

FragPipe's spectral library generation step: FragPipe's own `tools/speclib/gen_con_spec_lib.py` with the arguments of `CmdSpecLibGen`.

## Description

The task directory is FragPipe's group work directory.
It holds `psm.tsv`, `peptide.tsv`, the validated pepXMLs and `filelist_speclibgen.txt` (one LC-MS path per line).
`gen_con_spec_lib.py` converts the pepXMLs (`fragpipe-speclib convert --pepxml`, the FragPipe default) or, when the convert extras hold `--decoy_prefix`, the `psm.tsv` rows (`convertpsm`).
It then aligns retention times (automatic reference run, Biognosys/ciRT/Pierce iRT or a user file, with the ciRT retry), builds the EasyPQP library and writes `library.tsv` with `easypqp_lib_export('Spectronaut')`.
The process environment is FragPipe's: `PYTHONIOENCODING=utf-8`, `OMP_NUM_THREADS`, `FRAGPIPE_PARSE_MZML_CMD` and no `OPENMS_DATA_PATH`.

## Inputs

| Channel      | Type                              | Description                                                                                                |
| ------------ | --------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| meta         | val(map)                          | Group metadata                                                                                             |
| psm          | path (stageAs `psm.tsv`)          | The group's psm.tsv                                                                                        |
| peptide      | path (stageAs `peptide.tsv`)      | The group's peptide.tsv                                                                                    |
| pepxml_files | path                              | Validated `interact-*.pep.xml` or PTMProphet `*.mod.pep.xml`                                               |
| spectra      | path (stageAs `spectra/*`)        | LC-MS files of `filelist_speclibgen.txt`                                                                   |
| rt_cal_file  | path (stageAs `rt_calibration/*`) | User RT calibration file, `[]` otherwise                                                                   |
| im_cal_file  | path (stageAs `im_calibration/*`) | User IM calibration file, `[]` otherwise                                                                   |
| config_cli   | val(string)                       | `SPECLIB_RT_CAL`, `SPECLIB_IM_CAL`, `SPECLIB_CONVERT_ARGS`, `SPECLIB_LIBRARY_ARGS`, `SPECLIB_INTERMEDIATE` |
| fasta        | path                              | Protein database with decoys                                                                               |

## Outputs

| Channel             | Type             | Description                                                                  |
| ------------------- | ---------------- | ---------------------------------------------------------------------------- |
| library_tsv         | tuple(val, path) | `library.tsv`, the DIA-NN input                                              |
| easypqp_files       | tuple(val, path) | `easypqp_files/` with `cmds.txt` and alignment PDFs                          |
| filelists           | tuple(val, path) | `filelist_speclibgen.txt`, `filelist_easypqp_library.txt`                    |
| argv                | tuple(val, path) | `speclibgen_argv.txt`, the command one argument per line                     |
| intermediates       | tuple(val, path) | Calibration references, per-run PSM splits, EasyPQP intermediates (optional) |
| versions_speclibgen | topic: versions  | FragPipe-SpecLib version                                                     |
| versions_fragpipe   | topic: versions  | FragPipe version                                                             |

## Parameters

| Parameter              | Default                                           | Description                                                       |
| ---------------------- | ------------------------------------------------- | ----------------------------------------------------------------- |
| ext.args               | (none)                                            | Not accepted: gen_con_spec_lib.py takes positional arguments only |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory                                          |

## Container

The FragPipe image given with `--fragpipe_container` (its Python environment must have the EasyPQP dependency versions FragPipe-SpecLib was built with).

## Usage

Used in the `FRAGPIPE_SPECLIB` subworkflow (`subworkflows/local/fragpipe_speclib/main.nf`); `DIANNFRAGPIPE` quantifies DIA data with its `library.tsv`.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [SpecLib tutorial](https://fragpipe.nesvilab.org/docs/tutorial_speclib.html)
