# Z-DNA Genomic Analysis Pipeline

A comprehensive pipeline for analyzing Z-DNA regions in genomic data, including interval generation, perturbation testing, VCF annotation, mutation concentration calculation, and sequence generation from mutations.

## Features

- **Interval Generation**: Create and analyze genomic intervals
- **Perturbation Testing**: Statistical testing of genomic region overlaps
- **VCF Annotation**: Annotate VCF files with genomic features
- **Mutation Analysis**: Calculate mutation concentrations in regions
- **Sequence Generation**: Generate reference and mutant sequences with flanking regions

## Pipeline Components

### 1. Interval Generation
Generates random genomic intervals matching the size distribution of Z-DNA regions and counts intersections with genomic features.

**Key Functions**:
- `file_process()`: Processes interval lengths
- `read_coordinates()`: Reads and filters genomic coordinates
- `generate_and_count_intersections()`: Monte Carlo simulation of interval placement

### 2. Perturbation Testing
Performs statistical tests on overlaps between Z-DNA regions and other genomic features.

**Key Functions**:
- `fisher_method()`: Combines p-values using Fisher's method
- `p_value()`: Calculates empirical p-values
- `process_file()`: Processes super-enhancer files

### 3. VCF Annotation
Bash script for annotating VCF files with genomic regions using bedtools and tabix.

### 4. Mutation Concentration
Calculates mutation concentrations in genomic regions.

**Key Functions**:
- `clean_dataframe()`: Cleans and validates mutation data
- `calculate_mutation_concentration()`: Computes SNP concentrations

### 5. Sequence Generation
Generates reference and mutant sequences with flanking regions from mutation data.

**Key Functions**:
- `generate_mutant_sequences()`: Main sequence generation function
- `extract_reference_sequence()`: Gets reference sequence
- `extract_sequence_with_flanks()`: Gets sequence with flanking regions

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

### Input Files
- **Genomic Data**: BED files with genomic coordinates
- **Mutation Data**: Excel/TSV files with mutation information
- **Reference Genome**: FASTA file (hg38)
- **Motif Data**: Excel files with HOMER and MEME motifs

### Output
- Generated intervals (text files)
- Statistical results (Excel files)
- Annotated VCF files
- Mutation concentration calculations
- Sequence FASTA files (reference and mutant)

## Installation

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
