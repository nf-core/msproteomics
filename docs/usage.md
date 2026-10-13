# nf-core/msproteomics: Usage

## Introduction

nf-core/msproteomics is a mass spectrometry proteomics preprocessing pipeline supporting DIA, DDA LFQ, TMT label check, and generic FragPipe workflows.
This document describes how to prepare input data, configure the pipeline, and run each workflow type.

> **New to the pipeline?** Start with the [Quick Start Guide](quickstart.md) for a 5-minute introduction.
> For database preparation, see the [Database Preparation Guide](database_guide.md).
> For pipeline design details, see the [Workflow Architecture](architecture.md).

## Input Format

### CSV Samplesheet

The pipeline accepts a CSV samplesheet describing your samples and their raw data file locations.

```bash
--input samplesheet.csv
```

The samplesheet must contain the following columns:

| Column         | Description                                                                                                 | Required |
| -------------- | ----------------------------------------------------------------------------------------------------------- | -------- |
| `sample`       | Sample identifier; for TMT rows, the channel's sample name (`NA` excludes the channel)                      | Yes      |
| `spectra`      | Path or URI to a .raw, .mzML or .dia file or a Bruker .d directory, optionally .gz/.tar/.tar.gz/.zip        | Yes      |
| `condition`    | Experimental condition; fills the TMT annotation and MSstats condition fields (defaults to the sample name) | No       |
| `label`        | Isobaric label channel (e.g., TMT126, TMT127N); leave empty for LFQ/DIA                                     | No       |
| `fraction`     | Fraction number for fractionated experiments; leave empty for single-shot                                   | No       |
| `experiment`   | FragPipe experiment name (the "Experiment" column of the FragPipe manifest)                                 | No       |
| `bioreplicate` | FragPipe bioreplicate, an integer (the "Bioreplicate" column of the FragPipe manifest)                      | No       |
| `data_type`    | FragPipe data type: `DDA`, `DDA+`, `DIA`, `DIA-Quant`, `DIA-Lib` or `GPF-DIA` (the manifest's 4th column)   | No       |

Compressed spectra and `.dia` files are for `--mode diann` only.
`--mode fragpipe` takes `.raw`, `.mzML` and Bruker `.d` LC-MS files and stops with an error on anything else (FragPipe 24.0 itself reads no compressed files).
FragPipe workflows group files like FragPipe's `InputLcmsFile.getGroup()`.
When the workflow runs SpecLibGen (`speclibgen.run-speclibgen=true`) or DIA-NN (`diann.run-dia-nn=true`), every file is in one group, whatever its `experiment` and `bioreplicate`.
With `experiment` and `bioreplicate` set, the group is `<experiment>_<bioreplicate>`.
With only `bioreplicate` set, the group is `exp_<bioreplicate>`.
With only `experiment` set, the group is `<experiment>`.
With both blank, every file is in one group, which is the FragPipe default.
FragPipe replaces every character of `experiment` other than letters, digits and `_` with `_`.
Each group gets one Philosopher filter/report, one results directory and one IonQuant/TMT-Integrator experiment column.
The `condition` column does not affect grouping.
A blank `data_type` is guessed from the file path as FragPipe does: `DDA` if the path contains `dda` (any case), `_Q1.`, `_Q2.`, `_Q3.` or `_diatracer.`; else `DIA` if it contains `DIA`; else `DDA+` if it contains `WWA`; else `DDA`.
The data type selects the MSFragger parameter set (`fragger.params`, `fragger_dia.params` for DIA and DIA-Lib, `fragger_gpfdia.params`, `fragger_dda_plus.params`), the PercolatorOutputToPepXML mode (DDA only for `DDA`) and which files DIA-Umpire (not `DDA`, `DDA+` or `DIA-Quant`) and diaTracer (`DIA`, `DIA-Lib`, `GPF-DIA` `.d` files) process; their outputs are searched as `DDA`.
`DIA-Quant` files are not searched and do not go into the spectral library; they are only for DIA-NN quantification.

#### Spectral library and DIA-NN quantification (FragPipe pipeline mode)

With `speclibgen.run-speclibgen=true`, the pipeline runs FragPipe's own `tools/speclib/gen_con_spec_lib.py` with the arguments FragPipe builds (`CmdSpecLibGen`).
Its work directory holds `psm.tsv`, `peptide.tsv`, the validated pepXMLs (`interact-*.pep.xml`, or PTMProphet's `*.mod.pep.xml`) and `filelist_speclibgen.txt`, as FragPipe's does.
`speclibgen.convert-pepxml=true` (the FragPipe default) builds the library from the pepXMLs.
`speclibgen.convert-psm=true` builds it from `psm.tsv` rows and adds `--decoy_prefix <database.decoy-tag>`, the glyco option `--labile_mods` (`speclibgen.easypqp.labile_mode`) and `--max_glycan_q`.
Glyco modes need `convert-psm=true`; otherwise the run stops with FragPipe's message.
RT calibration (`speclibgen.easypqp.rt-cal`) is automatic, `Biognosys_iRT`, `ciRT`, `Pierce_iRT`, or `User provided RT calibration file` with the path in `speclibgen.easypqp.select-file.text`.
IM calibration is automatic or `User provided IM calibration file` with `speclibgen.easypqp.select-im-file.text`.
User calibration files must exist; they are staged into the task.

With `diann.run-dia-nn=true` and at least one `DIA` or `DIA-Quant` file, DIA-NN quantifies those files after SpecLibGen.
The pipeline runs the DIA-NN binary FragPipe 24.0 ships (`tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8`, DIA-NN 1.8.2 Beta 8) from `--fragpipe_container` and stops when the binary reports another version.
DIA-NN reads SpecLibGen's `library.tsv`, or `diann.library` when the workflow sets an external library.
The command is FragPipe's: `--lib`, `--threads`, `--verbose 1`, `--out dia-quant-output/report.tsv`, `--qvalue`, `--matrix-qvalue`, `--matrices`, protein inference (`--no-prot-inf`, or `--relaxed-prot-inf --fasta` with `diann.redo-protein-inference`), `--smart-profiling`, the quantification strategy, `--report-lib-info` for MSstats, plexDIA flags, `diann.cmd-opts`, and `--cfg filelist_diann.txt` with the inputs largest first.
A `DIA` `.raw` file is quantified from its mzML (MSFragger's `_uncalibrated.mzML` in FragPipe); `DIA-Quant` files are used as given.
After DIA-NN, FragPipe's `Propagation`, `DiannToMsstats` (when `diann.generate-msstats=true`) and `FragReporter` run, and `experiment_annotation.tsv` is written.
FragPipe keeps the blank group in its work-directory root, and so does this pipeline: its results are published at the top of each tool's output directory (for example `philosopher/psm.tsv`).
With every group blank, IonQuant runs without `--multidir`, and FragPipe refuses match-between-runs (`ionquant.mbr=1`) in that case; the pipeline stops with FragPipe's message.

Example samplesheet (DIA or DDA LFQ, one group per biological replicate):

```csv
sample,spectra,condition,label,fraction,experiment,bioreplicate
Sample_A1,/data/raw/A1.raw,control,,,control,1
Sample_A2,/data/raw/A2.raw,control,,,control,2
Sample_B1,/data/raw/B1.raw,treated,,,treated,1
Sample_B2,/data/raw/B2.raw,treated,,,treated,2
```

For TMT, one plex is one FragPipe group, so give every file of a plex the same `experiment` (and the same `bioreplicate`, or leave it blank).
List one row per channel and file; each file is converted and searched once.
For TMT rows, `sample` is the channel's sample name in the TMT annotation (FragPipe's `<channel> <sample>` annotation file), and `NA` excludes a channel.
Every channel of the label type set in the workflow file (`tmtintegrator.channel_num`) must be listed for each plex, as FragPipe requires.
A channel may appear once per plex; repeating it for each fraction file with the same sample and condition is allowed.
Sample names other than `NA` must be unique across plexes.
A plex without any labelled rows gets FragPipe's default annotation, which names each channel `<group>_<channel>` (`_<channel>` for the blank group).
No annotation is built when TMT-Integrator is disabled in the workflow file.

Example samplesheet (TMT-6, one plex split over two fractions):

```csv
sample,spectra,condition,label,fraction,experiment,bioreplicate
Bridge,/data/raw/pool1_f1.raw,reference,TMT126,1,Pool1,
Ctrl_1,/data/raw/pool1_f1.raw,control,TMT127N,1,Pool1,
Ctrl_2,/data/raw/pool1_f1.raw,control,TMT128C,1,Pool1,
Trt_1,/data/raw/pool1_f1.raw,treated,TMT129N,1,Pool1,
Trt_2,/data/raw/pool1_f1.raw,treated,TMT130C,1,Pool1,
NA,/data/raw/pool1_f1.raw,NA,TMT131N,1,Pool1,
Bridge,/data/raw/pool1_f2.raw,reference,TMT126,2,Pool1,
```

The remaining rows repeat the channels for `pool1_f2.raw` (see `assets/samplesheet_tmt.csv`).

The pipeline internally generates an SDRF file from the samplesheet via the `GENERATESDRFFROMSAMPLESHEET` module.
This SDRF is used for bookkeeping and compatibility with downstream tools that expect SDRF-formatted metadata.

## Analysis Modes

The `--mode` parameter selects the analysis engine:

| Mode       | Description                                 | Engine                                       |
| ---------- | ------------------------------------------- | -------------------------------------------- |
| `diann`    | DIA quantitative proteomics                 | [DIA-NN](https://github.com/vdemichev/DiaNN) |
| `fragpipe` | DDA LFQ, TMT, or generic FragPipe workflows | [FragPipe](https://fragpipe.nesvilab.org/)   |

### DIA Mode (`--mode diann`)

DIA workflows use DIA-NN for library-free or library-based DIA analysis.
DIA-specific variants (e.g., phospho) are applied via `-c conf/variants/*.config` files (see [DIA Method Variants](#dia-method-variants)).

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

`--database` or `--proteome` ([iProteomes](#iproteomes)) is required. In FragPipe mode the FASTA is searched as given and must already contain decoys tagged with the workflow `database.decoy-tag`; `--add_philosopher_decoy` / `--add_philosopher_contaminants` adds decoys and contaminants first (FragPipe's "Add decoys", `philosopher database --custom <fasta> [--nodecoys] [--contam --contamprefix]`).
To use a custom database, provide `--database /path/to/database.fasta`.

### FragPipe Mode (`--mode fragpipe`)

FragPipe mode supports DDA LFQ, TMT label check, TMT quantification, and generic FragPipe workflows.
The specific workflow is controlled by the `--fragpipe_workflow` and `--tmt_mode` parameters.

**DDA LFQ** -- uses MSFragger, MSBooster, Percolator, ProteinProphet, Philosopher, and IonQuant:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_container docker.io/<your-org>/fragpipe:24.0 \
  --fragpipe_workflow LFQ-MBR.workflow \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**TMT Label Check** -- assess TMT labeling efficiency:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_container docker.io/<your-org>/fragpipe:24.0 \
  --tmt_mode labelcheck \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**TMT Quantification** -- full TMT isobaric quantification:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_container docker.io/<your-org>/fragpipe:24.0 \
  --tmt_mode quant \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**Generic FragPipe** -- run any FragPipe workflow by providing a `.workflow` configuration file:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_container docker.io/<your-org>/fragpipe:24.0 \
  --fragpipe_workflow /path/to/custom.workflow \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

### Understanding .workflow Files

FragPipe `.workflow` files are key=value configuration files that control which tools run and their parameters.
The pipeline ships three built-in workflow files in `assets/`:

| File                          | Use Case                                              | Key Tools                                                  |
| ----------------------------- | ----------------------------------------------------- | ---------------------------------------------------------- |
| `LFQ-MBR.workflow`            | DDA label-free quantification with match-between-runs | MSFragger, MSBooster, Percolator, ProteinProphet, IonQuant |
| `TMT-labelcheck-229.workflow` | TMT label check for TMT6/10/11 (mass 229 Da)          | MSFragger, Percolator, ProteinProphet                      |
| `TMT-labelcheck-304.workflow` | TMT label check for TMTpro/16/18 (mass 304 Da)        | MSFragger, Percolator, ProteinProphet                      |

To create a custom `.workflow` file:

1. Open the [FragPipe GUI](https://fragpipe.nesvilab.org/) and configure your analysis
2. Save the workflow configuration (File > Export Workflow)
3. Provide the exported file to the pipeline with `--fragpipe_workflow /path/to/custom.workflow`

Key parameters users commonly tune in `.workflow` files:

- **Enzyme**: trypsin, lysc, chymotrypsin, etc.
- **Missed cleavages**: typically 1-2
- **Variable modifications**: oxidation (M), acetylation (protein N-term)
- **Mass tolerances**: precursor and fragment mass tolerances
- **Calibrate mass**: 0 (disabled), 1 (single pass), 2 (two pass, recommended)

> [!NOTE]
> TMT label check workflows use TMT as a **variable** modification (not fixed) to detect unlabeled peptides.
> This is intentional and differs from standard TMT quantification workflows.

FragPipe searches the FASTA the workflow names in `database.db-path`; the pipeline searches `--database` (or the iProteome), which replaces the workflow's `database.db-path`.
The workflow's path is a path on the machine that saved it (it cannot be read from Tower or AWS Batch), and renamed copies of one FASTA are normal, so a different file name only logs a warning naming both files.
The run stops, as FragPipe does ("FASTA file path is empty or the file is corrupted"), when `database.db-path` is FragPipe's `databases` folder, the value of a second-pass workflow, and on an empty `database.decoy-tag` ("Decoy protein prefix is empty. Please fix it in Database tab.").

Run-tab settings (`tab-run.*`) in pipeline mode:

- `tab-run.export_matched_fragments`: runs the matched-fragment export, as FragPipe.
- `tab-run.write_sub_mzml=true` runs FragPipe's first pass of its two-pass search, as FragPipe: MSFragger writes `<run>.mzBIN_all` and gets `--no-ptm-score-in-pin`, Philosopher Report prints decoys (the Run tab turns "Print decoys" on), Percolator saves `<run>_percolator.weights`, WriteSubMzml writes `fragpipe/<run>_sub.mzML`, and `fragpipe/fragpipe-second-pass.workflow` and `fragpipe/fragpipe-files-second-pass.fp-manifest` are written for the second pass.
  As in FragPipe, the run stops when MSFragger is off or a manifest file is DIA, GPF-DIA, DIA-Lib or DDA+; it also stops with a split database search (`msfragger.misc.slice-db` > 1), where FragPipe hands `--no-ptm-score-in-pin` to `msfragger_pep_split.py` as an LC-MS file and MSFragger stops with `Could not find MS file at <workdir>/--no-ptm-score-in-pin` (FragPipe 24.0 headless run).
  The second pass runs as in FragPipe on a samplesheet of the `_sub.mzML` files with `fragpipe-second-pass.workflow`: for a `_sub.mzML` input with its `<run>_percolator.weights` anywhere below the input's folder, MSFragger gets `--no-ptm-score-in-pin` and Percolator `--init-weights <weights> --static --override`; two weights files of the same name stop the run.
  As FragPipe, the second-pass workflow has an empty `database.decoy-tag` and the FragPipe `databases` folder as database, so the run stops on the empty decoy tag until both are set.
  MSFragger runs once per LC-MS file in the pipeline, so FragPipe's split of one MSFragger command at the command-line length limit (several commands per data type, each deciding `--no-ptm-score-in-pin` for its own files) does not arise.
- `tab-run.sub_mzml_prob_threshold`: WriteSubMzml leaves out the scans whose PSM probability in `psm.tsv` is above it.
- `tab-run.delete_temp_files` is ignored. FragPipe deletes temporary files from its work folder after the run (FragpipeRun.java:547-578; its own message says "It won't affect the result"); the pipeline keeps intermediate files in the Nextflow work directory and publishes the same results either way.

Settings after which FragPipe 24.0 aborts the run stop the pipeline with FragPipe's message:

- MSBooster with an open search (Da units, lower < -3, upper > 3), with neither Predict RT nor Predict spectra, or with a non-DIA-NN model (or Find best model) and an empty `msbooster.koina-url` (`CmdMSBooster.java:106-131`).
- Percolator without pin and pepXML, or PeptideProphet without pepXML, in `msfragger.output_format` (`FragpipeRun.java:252-266`).
- A labile search mode with `remainder_fragment_masses` and `localize_delta_mass=false` (`CmdMsfragger.java:418-427`).
- An enabled variable modification whose sites do not match `[A-Znc\[\]\^\*]+` (`CmdMsfragger.java:331-335`).
- A malformed `ptmshepherd.varmod_masses`, also with PTM-Shepherd off (`PtmshepherdPanel.java:816-829`).
- An IonQuant `excludemods` or `formula` value with whitespace, or a `light`, `medium` or `heavy` value with whitespace unless `ionquant.use-labeling=false` (FragPipe skips the labels then, `CmdIonquant.java:296-301`): FragPipe passes each as one argument, the pipeline splits the IonQuant options on spaces (the pipeline's own message, not FragPipe's).

Files named by a workflow (`msbooster.spectral-library-path`, `opair.oxonium_filtering_file`, `ptmshepherd.annotation_file` with the custom or glyco annotation option, the speclibgen RT and IM calibration files) are staged into the task; a path that does not exist stops the run.
`diann.cmd-opts` is passed to DIA-NN as one argument, as FragPipe does.

The `--fragpipe_mode` parameter controls execution mode:

- `pipeline` (default): Modular Nextflow subworkflow execution
- `headless`: Single-process FragPipe headless execution

### SAINTexpress

SAINTexpress is an optional post-run step, not part of FragPipe's run: FragPipe runs it only from its Downstream tab's own Run button (`TabDownstream.java:114-145,220`), and `FragpipeRun.java` has no SAINTexpress step.
It runs when the `.workflow` file sets `saintexpress.run-saint-express=true`; there is no pipeline flag.
The pipeline runs it after the main analysis, as the Downstream tab's Run button does, logs that it does so, and writes it to its own output directory, `downstream/saintexpress/`: both modes (`spc`, `int`) on the work dir's `reprint.<mode>.tsv`.
IonQuant writes those tables when it runs.
Otherwise the pipeline runs Philosopher Abacus with `--reprint --razor` to produce them, which FragPipe's own run conditions would not do (user decision, parity review Appendix C); the Abacus outputs then appear only in the pipeline results.
A mode whose table is missing is skipped with FragPipe's error `Input file <path> is not accessible.`, logged and written to `downstream/saintexpress/saintexpress-<mode>/skipped_mode.txt`.
SAINTexpress needs control IPs to score well, but FragPipe does not require them, and neither does the pipeline.
A bait is a control when its name starts with `CONTROL` (`CmdSaintExpress.java:117-123`); Philosopher's reprint table labels a group `CONTROL` when its name contains `control` in any case (`lib/aba/pro.go:102-122`).
With fewer controls than `saintexpress.virtual-controls` (`-L`, 100 by default) SAINTexpress prints `L is larger than the number of control IPs.` and continues.
With no controls and only one or two IPs it can abort (exit 134, `terminate called after throwing an instance of 'char const*'` or `nlopt::roundoff_limited`): FragPipe 24.0's `SAINTexpress-int` and `SAINTexpress-spc` do this on such input, and FragPipe's Downstream Run fails the same way.
The task then fails with SAINTexpress's output and names the step; add control experiments (names starting with `CONTROL`) or set `saintexpress.run-saint-express=false`.

## Key Parameters

| Parameter             | Description                                                | Default                                  |
| --------------------- | ---------------------------------------------------------- | ---------------------------------------- |
| `--input`             | Path to CSV samplesheet                                    | Required                                 |
| `--mode`              | Analysis mode: `diann` or `fragpipe`                       | Required                                 |
| `--outdir`            | Output directory                                           | Required                                 |
| `--database`          | FASTA protein database                                     | `--database` or `--proteome` is required |
| `--proteome`          | iProteome key ([iProteomes](#iproteomes))                  | None                                     |
| `--fragpipe_workflow` | FragPipe `.workflow` file or workflow name (FragPipe mode) | `assets/LFQ-MBR.workflow`                |
| `--tmt_mode`          | TMT analysis mode: `labelcheck` or `quant` (FragPipe mode) | None                                     |
| `--fragpipe_mode`     | FragPipe execution mode: `pipeline` or `headless`          | `pipeline`                               |

## Container Requirements

### DIA Workflows

DIA workflows use publicly available containers for DIA-NN and associated tools.
No special licensing is required.

### FragPipe-Based Workflows (DDA LFQ, TMT, Generic)

Every FragPipe tool runs from one licensed FragPipe image, given with `--fragpipe_container` (there is no default image), as FragPipe runs all tools from one distribution.
With a container engine and no `--fragpipe_container`, the run stops.
You must build that image with the licensed tools (MSFragger, IonQuant, diaTracer) under your own license terms.
See [Docker build instructions](fragpipe-docker/README.md) for academic and commercial container builds.

**FragPipe is free for academic use.**
Commercial users must obtain a license from Fragmatics.

## Profiles

Use `-profile` to select a software packaging method:

| Profile        | Description                                    |
| -------------- | ---------------------------------------------- |
| `docker`       | Run with Docker containers (recommended)       |
| `singularity`  | Run with Singularity containers                |
| `apptainer`    | Run with Apptainer containers                  |
| `podman`       | Run with Podman containers                     |
| `conda`        | Run with Conda environments (last resort)      |
| `test`         | Minimal stub test for CI validation (DIA mode) |
| `test_dia`     | DIA workflow test with small DIA dataset       |
| `test_dda_lfq` | DDA LFQ workflow test with small DDA dataset   |
| `test_tmt`     | TMT label check workflow test                  |
| `test_tmtq`    | TMT quantification workflow test               |

Multiple profiles can be combined: `-profile test_dia,docker`

## Running the Pipeline

### Basic Execution

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

### Using a Parameters File

```bash
nextflow run nf-core/msproteomics -profile docker -params-file params.yaml
```

```yaml
input: "samplesheet.csv"
mode: "fragpipe"
fragpipe_container: "docker.io/<your-org>/fragpipe:24.0"
fragpipe_workflow: "LFQ-MBR.workflow"
database: "/path/to/uniprot_human.fasta"
outdir: "results"
```

### Resuming a Failed Run

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker \
  -resume
```

### Updating the Pipeline

```bash
nextflow pull nf-core/msproteomics
```

### Reproducibility

Specify a pipeline version with `-r` to ensure reproducible results:

```bash
nextflow run nf-core/msproteomics -r 1.0.0 \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

## Configuration Hierarchy

The pipeline loads configuration files automatically based on the `--mode` parameter:

**DIA workflows** (`--mode diann`):

1. `nextflow.config` -- root config, always loaded
2. `conf/base_configs/diann.config` -- DIA-NN computed params and helper functions
3. `conf/variants/*.config` -- variant-specific overrides (e.g., diaphos) applied via `-c`
4. `conf/instruments/*.config` -- instrument-specific overrides applied via `-c`

**FragPipe workflows** (`--mode fragpipe`):

1. `nextflow.config` -- root config, always loaded (all params defined here)
2. `--tmt_type` selects the TMT plex (TMT6, TMT10, TMT11, TMT16, TMT18, TMTPRO); the correct labelcheck workflow file is auto-selected
3. `conf/instruments/*.config` -- instrument-specific overrides applied via `-c`

## Custom Configuration

### Resource Requests

To change compute resources for specific processes, see the [nf-core resource tuning documentation](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources).

### Custom Tool Arguments

To pass additional arguments to specific tools, use the `ext.args` directive in a custom config file.
See the [nf-core tool arguments documentation](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments).

### Resource Recommendations

Resource requirements depend on dataset size and analysis type:

| Dataset Size | Files  | Recommended Memory | Notes                                  |
| ------------ | ------ | ------------------ | -------------------------------------- |
| Small        | < 20   | 16 GB              | Default resources sufficient           |
| Medium       | 20-100 | 32 GB              | Increase MSFragger and IonQuant memory |
| Large        | 100+   | 64 GB              | Consider splitting into batches        |

- **DIA-NN library generation** is the most memory-intensive step; allocate 32-64 GB for large proteomes or phosphoproteomics
- **MSFragger** memory scales with database size; large databases (>100,000 sequences) may need 32+ GB
- **FragPipe headless mode** runs all tools in a single process and requires memory sufficient for the largest tool
- **FragPipe pipeline mode** distributes work across processes, allowing smaller per-process memory allocation

For cloud environments, memory-optimized instances are recommended (e.g., AWS r6i family, GCP n2-highmem).

> [!WARNING]
> Do not use `-c <file>` to specify pipeline parameters.
> Custom config files specified with `-c` must only be used for resource tuning, output directories, or module arguments (`ext.args`).

### Institutional Configs

The pipeline dynamically loads configurations from [nf-core/configs](https://github.com/nf-core/configs).
Check if your institution has a pre-configured profile available.

## Running on Seqera Platform (Tower)

### Nextflow Version

Set the Nextflow version in your launch pre-run script:

```bash
export NXF_VER=25.10.4
```

### Compute Environments

FragPipe processes require compute environments with sufficient memory.
Recommended instance families: `r6i` (memory-optimized).

## Pre-built DIA-NN Spectral Library

In-silico library generation (`DIANN_INSILICOLIBRARYGENERATION`) is the slowest DIA-NN step.
Give an existing library with `--diann_speclib` to skip it together with `QUANTMSUTILS_DIANNCFG`:

```bash
nextflow run main.nf --mode diann --input samplesheet.csv --database database.fasta \
    --diann_speclib /path/to/library.speclib --outdir results -profile docker
```

- The library goes to DIA-NN with `--lib` for the preliminary analysis; the empirical library is still assembled from your runs.
- DIA-NN 1.8.1 (the pipeline image) reads the compact `.speclib` format and table libraries (`.tsv`, `.xls`, `.txt`, `.csv`). `.parquet` libraries need DIA-NN 1.9 or later.
- `--database` is still required: DIA-NN uses the FASTA for protein inference.
- To skip the preliminary analysis and empirical library assembly as well, set `--diann_skip_preliminary_analysis` with `--diann_empirical_assembly_lib` and `--diann_empirical_assembly_log` (the library and log of a previous assembly).

## iProteomes

iProteomes are versioned reference proteomes with hosted indices, modelled on nf-core iGenomes.
`--proteome <key>` selects an entry of `conf/iproteomes.config` (`human_uniprot_reviewed_2026_03`, `mouse_uniprot_reviewed_2026_03`, `yeast_uniprot_reviewed_2026_03`).
`--iproteomes_base` is the root directory of the store; the repository has no default for it.

```bash
nextflow run main.nf --mode diann --input samplesheet.csv \
    --proteome human_uniprot_reviewed_2026_03 \
    --iproteomes_base /path/to/iproteomes --outdir results -profile docker
```

Entry paths in `conf/iproteomes.config` are relative to `--iproteomes_base` and joined with it at run time, so any params source (command line, `-params-file`, `-c`) works.

Each entry has:

| Field            | Meaning                                                                                       |
| ---------------- | --------------------------------------------------------------------------------------------- |
| `organism`       | Organism of the SDRF bookkeeping when `--organism` is not set                                 |
| `fasta`          | Protein FASTA                                                                                 |
| `fragpipe_fasta` | Optional FASTA with `rev_` decoys and contaminants (`philosopher database --custom --contam`) |
| `diann_speclib`  | Optional DIA-NN predicted spectral library of `fasta`                                         |

Resolution follows iGenomes with `--star_index`: a user parameter wins over the iProteome, which wins over building in the pipeline.

| Mode     | User parameter                  | iProteome                                             | Built in the pipeline                                                                                                                                                   |
| -------- | ------------------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| DIA-NN   | `--database`, `--diann_speclib` | `fasta`; `diann_speclib` when `--database` is not set | In-silico library (`DIANN_INSILICOLIBRARYGENERATION`)                                                                                                                   |
| FragPipe | `--database`                    | `fragpipe_fasta`, else `fasta`                        | Decoys and contaminants only with `--add_philosopher_decoy` / `--add_philosopher_contaminants`; otherwise `CHECKDBCONFIG` fails like FragPipe on a FASTA without decoys |

- `--database` replaces the iProteome FASTA, so the hosted DIA-NN library (built from the iProteome FASTA) is not used.
- The hosted DIA-NN library is used only when the run's digest and modification settings (`--enzyme`, `--fixed_mods`, `--variable_mods`, `--max_mods`, `--allowed_missed_cleavages`, peptide length, precursor charge, `--met_excision`, the m/z ranges, the library mass accuracies and `--additional_insilico_args`) equal those recorded in the iProteome's `<iproteomes_base>/<key>/meta.json` (`diann_library_settings`, derived from `build_commands.diann`; numbers compare as numbers, modification lists as sets); otherwise the run stops and names the differing settings.
- With `--add_philosopher_decoy` / `--add_philosopher_contaminants` the plain `fasta` is updated by Philosopher and `fragpipe_fasta` is not used.
- An entry that declares an index (`diann_speclib`, `fragpipe_fasta`) whose path is empty, missing or unreadable stops the run; an entry that does not declare it has the pipeline build it.
- `--save_reference` publishes what the pipeline built to `<outdir>/proteome/`: the DIA-NN library to `diann/`, the Philosopher-updated database to `fragpipe/`.
- `--iproteomes_ignore` skips `conf/iproteomes.config`.
- To add an entry, copy one in `conf/iproteomes.config` (or a custom `-c` file) with paths relative to `--iproteomes_base`, and put a `meta.json` in `<iproteomes_base>/<key>/` for a hosted DIA-NN library.

## DIA-NN Mass Accuracy

With `--diann_mass_acc_automatic true` (default) DIA-NN takes the mass accuracies its preliminary analysis recommends.
With `--diann_mass_acc_automatic false` it uses `--precursor_mass_tolerance` (MS1, default 5 ppm) and `--fragment_mass_tolerance` (MS2, default 18 ppm); DIA-NN supports ppm only.
`--diann_library_ms1_acc` / `--diann_library_mass_acc` set the accuracy of in-silico library generation.

## Instrument-Specific Settings

DIA-NN in-silico library generation uses different m/z ranges and mass accuracy settings depending on the mass spectrometer.
The pipeline ships per-instrument config files in `conf/instruments/` that override the defaults when included with `-c`.
Each also sets `precursor_mass_tolerance` and `fragment_mass_tolerance` to the library MS1 and MS2 accuracies, which DIA-NN uses when `--diann_mass_acc_automatic false`.

| Abbreviation | Instrument          | Config File                                       | diann_min_pr_mz | diann_max_pr_mz | diann_min_fr_mz | diann_max_fr_mz | diann_library_mass_acc | diann_library_ms1_acc |
| ------------ | ------------------- | ------------------------------------------------- | --------------- | --------------- | --------------- | --------------- | ---------------------- | --------------------- |
| ASC          | Thermo Ascend       | `conf/instruments/thermo_ascend.config` (default) | 350             | 1050            | 200             | 1800            | 18                     | 5                     |
| FLX          | Bruker timsTOF Flex | `conf/instruments/bruker_flex.config`             | 100             | 1700            | 100             | 1700            | 15                     | 15                    |

Usage example:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/instruments/bruker_flex.config \
  -profile docker
```

Thermo Ascend users can omit the `-c` flag since the defaults are already optimized for Ascend instruments.

## DIA Method Variants

Specialized DIA analysis modes are available as variant config files in `conf/variants/`.
These configs adjust DIA-NN parameters (e.g., variable modifications, phosphosite monitoring) for specific experimental designs.

| Variant               | Config File                    | Description                                        |
| --------------------- | ------------------------------ | -------------------------------------------------- |
| DIA Phosphoproteomics | `conf/variants/diaphos.config` | Phosphosite localization with UniMod:21 monitoring |

Usage example:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/variants/diaphos.config \
  -profile docker
```

Variant configs can be combined with instrument configs by specifying multiple `-c` flags:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/instruments/bruker_flex.config \
  -c conf/variants/diaphos.config \
  -profile docker
```

## Nextflow Memory Requirements

To limit Nextflow's JVM memory usage, add the following to your environment:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Golden FragPipe tests

Golden tests check that a module reproduces FragPipe 24.0 exactly: the module gets the exact input files FragPipe used, and its command line and every output file are compared with FragPipe's.
They need licensed-environment data (a captured FragPipe run and the licensed FragPipe image), so they run outside public CI with the nf-test tag `golden`.

### Capture a FragPipe run

Run FragPipe 24.0 headless, keeping every intermediate file, then capture it:

```bash
python3 bin/golden_capture.py \
    --log <run>/work/log_<timestamp>.txt \
    --dataset <dataset> \
    --workflow TMT-labelcheck-304_perfile \
    --root /path/to/golden
```

If the run dir holds empty placeholder spectra files (for example Docker bind-mount targets), pass `--lcms-source <dir with the real files>`; an empty spectra file without it is an error.
`--root` defaults to `$GOLDEN_ROOT` and is required when that is unset; an unwritable root is an error.
The capture writes `<root>/<dataset>/<workflow>/` with:

| Path                                       | Content                                                                                                                                                                                                                                                                   |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `manifest.json`                            | every executed FragPipe step: tool, normalized argv, work dir, exit code, parallel batch and concurrent steps, inputs (with `golden_state`), outputs (with final state: present, moved, deleted, overwritten), flags for anything inferred; sha256 of every captured file |
| `fragpipe_log.redacted.txt`                | the FragPipe log with license arguments and license lines removed                                                                                                                                                                                                         |
| `files/work/`, `files/input/`, `files/db/` | FragPipe work dir, spectra dir (symlinks resolved), database dir                                                                                                                                                                                                          |

The FragPipe log is never copied verbatim, license-named files are refused, and the capture aborts if a copied file references a license file or `--license` argument.
In the redacted log and the manifest, the value of any license flag (`--license`, `-licence`, `--License-key`, ..., with a space or `=`) becomes `{LICENSE}` in place and any `license.dat` path becomes `{LICENSE_PATH}`; the argument count is unchanged, command lines are never dropped, and only informational lines that mention a license are dropped.
Logs already scrubbed this way are accepted.

FragPipe runs the per-group steps of a parallel group concurrently, so their log output interleaves.
The manifest records each step's `parallel_batch` and `concurrent_with`, and flags steps whose `Process '<name>' finished` marker cannot be told apart from a concurrent step of the same name.

Some steps rewrite their inputs in place (IonQuant, PTM-Shepherd in glyco mode, MBG, `Rewrite pepxml`, `AppendFile`), so a captured run only holds the rewritten version.
Each input carries a `golden_state` (`available`, `fragpipe_install`, `rewritten_by_step:N`, `moved_by_step:N`, `deleted_by_step:N`), and `golden_inputs_available` is false when any input is gone.
Golden tests for such steps take their inputs from a run that stops just before the step, for example a run that ends after Philosopher Report.
`bin/fragpipe_log_to_manifest.py <log> -o manifest.json` builds the manifest alone.

### Run the golden tests

```bash
export FRAGPIPE_CONTAINER=<licensed FragPipe 24.0 image>   # the FragPipe tools of the golden tests run in it
export GOLDEN_ROOT=/path/to/golden                         # captured FragPipe runs; every golden test fails without it
export OPAIR_DET_DLL=/path/to/PTMLocalization.dll          # deterministic O-Pair 1.1.10 build (O-Pair golden test)
NXF_OFFLINE=true nf-test test --tag golden
```

Stock O-Pair 1.1.10 breaks ties between equally scored glycan placements in a random order, so its `psm.tsv` differs between identical runs.
The O-Pair golden test compares against `PXD017646/glyco-O-Pair_perfile_det`, a FragPipe run with a `PTMLocalization.dll` that builds its glycan boxes in enumeration order (the order of stock O-Pair with `-t 1`).
It mounts `$OPAIR_DET_DLL` over the image's `tools/opair/PTMLocalization.dll`, checks its md5 (`01fca5e12c419a45880c5b610bcd0175`) and requires `psm.tsv` to match exactly.

Each golden test has its own `golden.config` next to it: it runs the tools in `$FRAGPIPE_CONTAINER` with docker and pins the tested processes to the threads and memory of the golden FragPipe run.

Each test (`modules/**/tests/golden.nf.test`) selects its steps from the manifest, stages the golden inputs, runs the module, and calls `bin/golden_step_check.py`, which fails on any of:

- argv differences after path normalization (step work dir, FragPipe install and external inputs become placeholders; argument order and flags are compared as is);
- output files that differ under the comparison rules below, or that the module did not write;
- module files with no FragPipe counterpart, except module logs passed with `--ignore-actual`.

Outputs FragPipe deleted or overwrote later (for example Percolator's `*_psms.tsv`) have no golden copy and are reported as `UNVERIFIABLE`; a step with no output that could be compared fails (`NO_VERIFIED_OUTPUT`).
The module's command line is split with POSIX shell rules, so quoted and empty arguments are compared like FragPipe's.

### Comparison rules

`bin/golden_compare.py` applies exactly the rules measured in the FragPipe-vs-FragPipe determinism runs; every other byte must match.
The Philosopher tolerances are the measured envelope of its nondeterministic FDR tie handling (`bin/philosopher_tables.py`, determinism doc section "Philosopher filter/report nondeterminism").
There are two rule profiles, selected with `--profile` (default `dda`); `bin/golden_step_check.py` picks `dia` when the golden run's spectra are DIA.

Profile `dda` (TMT label check, DDA LFQ, glyco DDA):

| File                                                                                                     | Normalization                                                                                                                      |
| -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `log_*.txt`, anything under `?/`                                                                         | skipped                                                                                                                            |
| `*.pepindex`                                                                                             | first 8 bytes dropped, run paths rewritten                                                                                         |
| `*.prot.xml`                                                                                             | `<tmpdir>/<uuid>/`, run paths, ISO `date=`/`time=`                                                                                 |
| `*.pepXML`, `*.pep.xml`                                                                                  | run paths, ISO `date=`/`time=`                                                                                                     |
| `filter.log`                                                                                             | run paths, `time=` stripped, `Parsing` lines as a set, PSM/peptide/ion counts ±2, decoy ±1, threshold ±2e-4 (Philosopher envelope) |
| `fragger*.params`, `*.workflow`, `*.job`, `*.fp-manifest`, `filelist_*.txt`, `experiment_annotation.tsv` | run paths                                                                                                                          |
| `psm.tsv`                                                                                                | Philosopher envelope: `Qvalue` within 1.01/max(T_above,1), at most 1 boundary PSM between the two runs' `filter.log` thresholds    |
| `peptide.tsv`, `ion.tsv`, `protein.tsv`                                                                  | spectral counts may move by the boundary PSMs mapped to the entry                                                                  |
| everything else                                                                                          | none (raw bytes)                                                                                                                   |

Profile `dia` adds the rules measured on the DIA runs:

| File                                             | Normalization                                                                                                                                                                           |
| ------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MSBooster/spectraRT.predicted.bin`              | skipped                                                                                                                                                                                 |
| `*.pdf` (EasyPQP alignment plots)                | `/CreationDate (D:...)`                                                                                                                                                                 |
| DIA-NN `report.log.txt`                          | run paths, `Current date and time:` text, `[m:ss] ` prefixes                                                                                                                            |
| `MSBooster/spectraRT*.tsv`                       | rows sorted (header first)                                                                                                                                                              |
| `peptide.tsv`, `ion.tsv`                         | full Philosopher envelope: `Qvalue` rule; `ion.tsv` keyed on stripped sequence + `Charge` + `M/Z`, `Expectation`/`Prev AA`/`Next AA` dropped; `peptide.tsv` without `Prev AA`/`Next AA` |
| `msbooster_params.txt`, `easypqp_files/cmds.txt` | run paths                                                                                                                                                                               |

Compare two run directories, for example two FragPipe runs:

```bash
python3 bin/golden_compare.py dirs runA/work runB/work --profile dda \
    --expected-prefix runA --actual-prefix runB \
    --expected-tmp runA/jtmp --actual-tmp runB/jtmp
```

Missing and extra files are failures, and so is a comparison in which no file was compared.

### Coverage

```bash
python3 bin/golden_coverage.py <golden>/<dataset>/<workflow>/manifest.json
```

lists every FragPipe step of the run and the golden test that covers it.
Unit tests for the scripts: `python3 -m unittest discover -s tests/python`.

## Comparing against FragPipe

`bin/compare_fragpipe_results.py` measures how closely an msproteomics run reproduces a FragPipe GUI/headless run on the same data.
It needs only Python 3 (standard library).

```bash
python3 bin/compare_fragpipe_results.py \
    --test results/ \
    --reference /path/to/fragpipe_output/ \
    --outdir comparison/
```

Both directories are searched recursively, so the msproteomics layout (`philosopher/<exp>/`, `ionquant/`, `tmtintegrator/`) and the FragPipe layout (`<exp>/`, `combined_*.tsv`, `tmt-report/`) can be compared directly.
Per-experiment `psm.tsv`, `peptide.tsv`, `ion.tsv` and `protein.tsv` are matched by experiment directory name.
Use `--experiment-map TEST=REF ...` when experiment names differ, and `--exclude REGEX` when a file is found twice.
Combined tables (`combined_protein.tsv`, `combined_peptide.tsv`, `combined_modified_peptide.tsv`, `combined_ion.tsv`) and TMT-Integrator tables (`abundance_*.tsv`, `ratio_*.tsv`) are matched by file name.
A table present on only one side, or a missing required column, is an error (exit code 2) that names the file.

The tool reports, per experiment and for the combined tables:

| Metric         | Description                                                                                                                 |
| -------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Counts         | PSMs, peptides, modified peptides, ions and proteins (decoys excluded)                                                      |
| Overlap        | Jaccard, only-test and only-reference counts for peptide sequences, modified peptides, ions and protein IDs                 |
| Shared PSMs    | Percentage of shared spectra with the same peptide and modified peptide, and Probability/Hyperscore/Expectation correlation |
| Quantification | Pearson and Spearman of log2 intensities (raw and MaxLFQ), median log2 ratio and percentage within 2-fold                   |
| TMT            | Per-channel correlation of TMT-Integrator log2 values                                                                       |
| Parameters     | Differing keys between `fragger.params` and `fragpipe.workflow` (informational, not gated)                                  |

Modified peptides are keyed on `Assigned Modifications` (position, residue and mass rounded to `--mod-decimals`), because `Modified Peptide` notation differs between tables.
The verdict is PASS (exit code 0) only when every gated metric meets `--max-count-deviation` (default 0.01), `--min-jaccard` (default 0.98) and `--min-spearman` (default 0.99); otherwise it is FAIL (exit code 1).
Correlations with fewer than `--min-quant-n` shared values are not computed and are listed as skipped.
Results are written to `comparison.md` and `comparison.json` in `--outdir`.

Unit tests run with `python3 -m unittest discover -s tests/python -v`.

## Further Reading

- [Quick Start Guide](quickstart.md) -- Get running in 5 minutes
- [Output Documentation](output.md) -- Interpreting pipeline results
- [Database Preparation](database_guide.md) -- Choosing and preparing FASTA databases
- [Workflow Architecture](architecture.md) -- Pipeline design and workflow routing
- [Troubleshooting](troubleshooting.md) -- Common errors and solutions
