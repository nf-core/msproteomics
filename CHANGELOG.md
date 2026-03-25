# nf-core/msproteomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - Initial Release

### Added

- DIA workflow (DIA-NN) with MSstats and MaxLFQ quantification
- DDA LFQ workflow (FragPipe: MSFragger, MSBooster, Percolator, IonQuant)
- TMT Label Check workflow (FragPipe: TMT labeling efficiency QC)
- Generic FragPipe workflow (configurable via .workflow files)
- Multi-instrument support (ASC, FLX)
- SDRF input format support
- Pre-built database selection for human, mouse, and yeast
- nf-core native pipeline structure
