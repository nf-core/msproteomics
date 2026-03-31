# nf-core/msproteomics: Troubleshooting

## Samplesheet Errors

- **"ERROR ~ Validation of 'input' file failed"** — your samplesheet has incorrect column names.
  The required columns are: `sample`, `spectra`, `condition`, `label`, `fraction`.
- **"does not match" pattern errors** — sample names cannot contain spaces.
  The `spectra` column must point to files ending with `.raw`, `.mzML`, `.d`, or `.dia`.
- **Common mistake**: using FASTQ-style samplesheets from other nf-core pipelines.
  This pipeline expects mass spectrometry file paths, not FASTQ read pairs.

## Database Errors

- **"FASTA file not found"** — check that the `--database` path exists and is accessible from all compute nodes.
- **"No decoy sequences detected"** — FragPipe headless mode requires decoy sequences.
  The pipeline adds them automatically, but if you provide a pre-processed database, ensure `rev_` prefix decoys are included.
- **DIA-NN database selection** — the database is auto-selected for supported organisms (human, mouse, rat, yeast, zebrafish, drosophila, C. elegans, E. coli).
  For other organisms, `--database` is required.

## Container and Licensing Errors

- **"MSFragger not found" or "IonQuant not found"** — the public FragPipe image does not include licensed tools.
  You must build a custom container with your own license.
  See [Docker build instructions](fragpipe-docker/README.md).
- **"Permission denied" on container** — ensure the Docker socket is accessible (`docker info`) or use `-profile singularity`.
- **Singularity cache issues** — set `NXF_SINGULARITY_CACHEDIR` to a writable directory with sufficient space.

## Memory and Resource Errors

- **"java.lang.OutOfMemoryError: Java heap space"** — increase JVM memory for Nextflow itself:

  ```bash
  export NXF_OPTS='-Xms1g -Xmx4g'
  ```

- **"Process exceeded memory limit"** — increase process memory in a custom config:

  ```groovy
  process {
      withName: 'MSFRAGGER' {
          memory = '32.GB'
      }
  }
  ```

- DIA-NN library generation is the most memory-intensive step.
  For large proteomes, allocate 32-64 GB.
- See [nf-core resource tuning documentation](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) for general guidance.

## Resume and Caching

- **`-resume` does not re-use cached results** — ensure you run from the same directory with the same `work/` directory.
- Changed parameters invalidate the cache for affected processes.
- **"Unable to serialize" warnings during resume** — safe to ignore unless results are incorrect.

## DIA-NN Specific Issues

- **"No precursors identified"** — check that instrument m/z ranges match your data.
  Use `-c conf/instruments/bruker_flex.config` for timsTOF data.
- **Library generation produces empty library** — verify input mzML files are not empty and contain MS2 spectra.
- **Mass accuracy issues** — adjust `--diann_library_mass_acc` and `--diann_library_ms1_acc` for your instrument.

## FragPipe Specific Issues

- **Headless vs pipeline mode**: `--fragpipe_mode pipeline` (default) runs each tool as a separate Nextflow process for better parallelization and caching.
  The `headless` mode runs FragPipe as a single process, which is useful for debugging or when using complex `.workflow` files.
- **".workflow file not found"** — provide the full path to the `.workflow` file or use a built-in workflow name.
- **TMT annotation errors** — ensure your samplesheet has correct `label` column values matching TMT channel names (TMT126, TMT127N, TMT128C, etc.).
- **"Philosopher workspace" errors** — these indicate database preparation issues.
  Check that your FASTA file is valid and not corrupted.

## Platform-Specific Issues

- **HPC/Slurm**: set `executor = 'slurm'` in a custom config and ensure Singularity is available on compute nodes.
- **Cloud (AWS/GCP)**: use Seqera Platform for easier orchestration.
  Ensure instance types have sufficient memory (`r6i` family recommended for FragPipe).
- **Local execution**: monitor disk space in the `work/` directory.
  Large proteomics datasets can consume 100+ GB.

## Getting Help

- Check the [nf-core documentation](https://nf-co.re/docs/usage/) for general Nextflow/nf-core troubleshooting.
- Search or open an issue on [GitHub](https://github.com/nf-core/msproteomics/issues).
- Join the [nf-core Slack](https://nf-co.re/join/slack) and ask in the `#msproteomics` channel.
