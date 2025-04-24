# Z-DNA Genomic Analysis Pipeline 

A comprehensive pipeline for analyzing Z-DNA regions in genomic data, including interval generation, perturbation testing, VCF annotation, mutation concentration calculation, and sequence generation from mutations.

## Features

- **Interval Generation**: Create and analyze genomic intervals
- **Perturbation Testing**: Statistical testing of genomic region overlaps
- **VCF Annotation**: Annotate VCF files with genomic features
- **Mutation Analysis**: Calculate mutation concentrations in regions
- **Sequence Generation**: Generate reference and mutant sequences with flanking regions

## Requirements

### Python Packages
- pandas>=1.3.0
- numpy>=1.21.0
- scipy>=1.7.0
- tqdm>=4.60.0
- sortedcontainers>=2.4.0
- biopython (for sequence generation)

### System Tools
- bedtools
- tabix

### Conda Environment

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.8
  - pandas
  - numpy
  - scipy
  - tqdm
  - sortedcontainers
  - bedtools
  - tabix
  - biopython
