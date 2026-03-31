# nf-core/msproteomics: Documentation

nf-core/msproteomics is a Nextflow pipeline for mass spectrometry proteomics supporting DIA (DIA-NN), DDA label-free quantification, TMT isobaric labeling, and generic FragPipe workflows.

## Getting Started

- [Quick Start](quickstart.md) -- Get running in minutes with copy-paste commands
- [Usage Guide](usage.md) -- Complete parameter reference and configuration options
- [Output Documentation](output.md) -- Description and interpretation of pipeline outputs

## Guides

- [Database Preparation](database_guide.md) -- Choosing and preparing FASTA protein databases
- [Workflow Architecture](architecture.md) -- Pipeline design, workflow routing, and subworkflow composition
- [Troubleshooting](troubleshooting.md) -- Common errors and solutions

## Reference

- [Module Documentation](modules/README.md) -- Detailed documentation for all 49+ analysis modules
- [FragPipe Docker Build](fragpipe-docker/README.md) -- Building custom FragPipe containers with licensed tools
- [Citations](../CITATIONS.md) -- References for all tools and methods used

## Example Samplesheets

- [DIA / DDA LFQ samplesheet](../assets/samplesheet.csv) -- Basic samplesheet for DIA or label-free experiments
- [TMT samplesheet](../assets/samplesheet_tmt.csv) -- TMT experiment with labeled channels and fractions

## Additional Resources

- [nf-core website](https://nf-co.re) -- General nf-core documentation, installation, and configuration guides
- [Nextflow documentation](https://www.nextflow.io/docs/latest/) -- Nextflow language and runtime reference
- [FragPipe documentation](https://fragpipe.nesvilab.org/) -- FragPipe tools and workflow configuration
- [DIA-NN documentation](https://github.com/vdemichev/DiaNN) -- DIA-NN analysis engine
