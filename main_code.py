# requirements:

# pandas>=1.3.0
# numpy>=1.21.0
# scipy>=1.7.0
# tqdm>=4.60.0
# sortedcontainers>=2.4.0

# environment.yml for conda:
# yaml

# channels:
#   - conda-forge
#   - bioconda
# dependencies:
#   - python=3.8
#   - pandas
#   - numpy
#   - scipy
#   - tqdm
#   - sortedcontainers
#   - bedtools
#   - tabix

# Standard library imports
import os
import sys
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Third-party imports
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import chi2
from sortedcontainers import SortedList



"""
Z-DNA GENOMIC ANALYSIS PIPELINE
==============================

A complete, reproducible pipeline for analyzing Z-DNA regions with:
1. Interval generation
2. Perturbation testing
3. Region grouping
4. VCF annotation
5. Mutation concentration calculation

All code is thoroughly commented and organized for:
- Transparency
- Reproducibility
- Maintainability
"""

# ======================
# 1. INTERVALS GENERATION
# ======================

def file_process(df):
    """
    Process DataFrame to extract interval lengths between start and end positions.
    
    Parameters:
        df (pd.DataFrame): Input dataframe with 'start' and 'end' columns
        
    Returns:
        list: List of interval lengths (end - start)
    """
    intervals = []
    # Using tqdm for progress tracking
    for _ in tqdm(range(len(df)), desc="Processing intervals"):
        row = df.iloc[_]
        start = row['start']
        end = row['end']
        length = end - start
        intervals.append(length)
    return intervals

def read_coordinates(file_path, chr):
    """
    Read genomic coordinates from BED file and filter by chromosome.
    
    Parameters:
        file_path (str): Path to BED file
        chr (str): Target chromosome (e.g., 'chr1')
        
    Returns:
        pd.DataFrame: Filtered and sorted dataframe with:
            - chrom: chromosome name
            - start: start position
            - end: end position
            - length: interval length
    """
    # Read BED file with specific column names and types
    coordinates = pd.read_csv(
        file_path,
        sep='\t',
        lineterminator='\n',
        header=None,
        names=['chrom','start', 'end'],
        dtype={'start': int, 'end': int},
        usecols=[0, 1, 2]
    )
    
    # Calculate interval lengths
    coordinates['length'] = coordinates['end'] - coordinates['start']
    
    # Filter by chromosome and sort by length
    coordinates = coordinates.loc[(coordinates['chrom'] == chr)]
    coordinates.sort_values(by='length', inplace=True)
    
    return coordinates

def subtract_interval(available_intervals, new_interval):
    """
    Subtract new interval from available intervals using interval arithmetic.
    
    Parameters:
        available_intervals (list): List of (start, end) tuples
        new_interval (tuple): Interval to subtract (start, end)
        
    Returns:
        list: Updated intervals after subtraction
    """
    updated_intervals = []
    for start, end in available_intervals:
        # Check for no overlap
        if not (end <= new_interval[0] or start >= new_interval[1]):
            # Handle partial overlaps
            if start < new_interval[0]:
                updated_intervals.append((start, new_interval[0]))
            if end > new_interval[1]:
                updated_intervals.append((new_interval[1], end))
        else:
            # Keep intervals that don't overlap
            updated_intervals.append((start, end))
    return updated_intervals

def generate_and_count_intersections(intervals, coordinates, num_iterations):
    """
    Generate random intervals and count intersections with genomic features.
    
    Parameters:
        intervals (list): List of interval lengths to generate
        coordinates (pd.DataFrame): Genomic coordinates dataframe
        num_iterations (int): Number of Monte Carlo iterations
        
    Returns:
        tuple: (total_results, errors)
            total_results: List of generated intervals per iteration
            errors: List of intervals that couldn't be placed
    """
    total_results = []
    errors = []

    # Main iteration loop with progress bar
    for _ in tqdm(range(num_iterations), desc="Iterations"):
        generated_intervals = []
        available_intervals = [(row['start'], row['end']) for _, row in coordinates.iterrows()]

        # Process each target interval length
        for interval in tqdm(intervals, desc="Processing Intervals", leave=False):
            valid_coordinate_found = False
            random.shuffle(available_intervals)  # Randomize selection
            
            # Try to place interval in available genomic space
            for start, end in available_intervals:
                if interval <= end - start:
                    new_start = random.randint(start, end - interval)
                    new_end = new_start + interval
                    new_interval = (new_start, new_end)
                    
                    generated_intervals.append(new_interval)
                    available_intervals = subtract_interval(available_intervals, new_interval)
                    valid_coordinate_found = True
                    break
            
            if not valid_coordinate_found:
                errors.append(interval)

        total_results.append(generated_intervals)
        
    return total_results, errors

def main(bed_data_2, output, file_path_2, chr):
    """
    Main execution function for interval generation.
    
    Parameters:
        bed_data_2: Input BED data
        output: Output file path
        file_path_2: Genomic coordinates file
        chr: Target chromosome
    """
    # Initialize data
    intervals = file_process(bed_data_2)
    coordinates = read_coordinates(file_path_2, chr)
    
    # Parallel processing setup
    num_processes = os.cpu_count()  
    num_iterations = 1000
    iterations_per_process = num_iterations // num_processes
    print(f"Using {num_processes} processes with {iterations_per_process} iterations each.")

    # Parallel execution with ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [
            executor.submit(
                generate_and_count_intersections,
                intervals,
                coordinates,
                iterations_per_process
            ) for _ in range(num_processes)
        ]
        
        # Initialize output file
        with open(os.path.expanduser(output), 'w') as f:
            f.write("Results:\n")
        
        # Process completed futures
        results = []
        errors = []
        completed_futures = []

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
                result, error = future.result(timeout=None)
                completed_futures.append(future)
                with open(os.path.expanduser(output), 'a') as f:
                    f.write(f'Intervals = {result}\n')
                    f.write(f'Errors from process: {error}\n')

        # Handle any remaining futures
        for future in futures:
            if future not in completed_futures:
                try:
                    result, error = future.result(timeout=None)
                    results.append(result)
                    errors.append(error)
                    with open(os.path.expanduser(output), 'a') as f:
                        f.write(f'Intervals = {result}\n')
                        f.write(f'Errors from process: {error}\n')
                except Exception as e:
                    print(f"An error occurred: {e}")

    # Finalize output
    with open(os.path.expanduser(output), 'a') as f:
        f.write("All processes completed.\n")

# Chromosome-wise execution
for i in range(24):
    if i == 0:
        pass  # Skip chromosome 0
    else:
        number = str(i)
        chr = 'chr' + number
        bed_data_2 = read_coordinates('ZDNABERT_hg38_generations.bed', chr)
        file_path_2 = 'hg38_data.bed'
        output = 'whole_genome_generated_z_dna_on_' + chr + '.txt'
        if __name__ == "__main__":
            main(bed_data_2, output, file_path_2, chr)


# ======================
# 2. PERTURBATION TESTING
# ======================

def fisher_method(p_values):
    """
    Combine p-values using Fisher's method.
    
    Parameters:
        p_values (array-like): List of p-values to combine
        
    Returns:
        float: Combined p-value
    """
    k = len(p_values)
    chi2_statistic = -2 * np.sum(np.log(p_values))
    return chi2.sf(chi2_statistic, df=2 * k)

def is_intersect(interval1, interval2):
    """
    Check if two genomic intervals intersect with minimum 10bp overlap.
    
    Parameters:
        interval1, interval2: Tuples of (start, end)
        
    Returns:
        bool: True if intervals intersect sufficiently
    """
    # Check minimum length requirement
    if interval1[1] - interval1[0] < 10 or interval2[1] - interval2[0] < 10:
        return False 
    
    # Check for any overlap
    if interval1[1] <= interval2[0] or interval2[1] <= interval1[0]:
        return False  
    
    # Calculate overlap
    start_intersect = max(interval1[0], interval2[0])
    end_intersect = min(interval1[1], interval2[1]) 
    return end_intersect - start_intersect >= 10

def p_value(real_data, combined_intervals, se_file, chr):
    """
    Calculate empirical p-value for overlap between real and simulated data.
    
    Parameters:
        real_data: Observed intervals
        combined_intervals: Simulated intervals
        se_file: Super-enhancer file
        chr: Chromosome
        
    Returns:
        float: Empirical p-value
    """
    fixed_intervals = process_fixed_intervals(read_coordinates(se_file, chr))
    intersections_real = optimize_intersection_check(fixed_intervals, real_data)
    total_results = []

    # Compare against all simulated intervals
    for iteration in tqdm(combined_intervals):
        list_1 = fixed_intervals
        list_2 = [node for node in iteration]
        intersections_count = optimize_intersection_check(list_1, list_2)
        total_results.append(intersections_count)            

    # Calculate empirical p-value
    total = sum(1 for i in total_results if i >= intersections_real)
    return total/len(combined_intervals)

def parse_intervals(file_path):
    """
    Parse interval data from output files.
    
    Parameters:
        file_path: Path to interval file
        
    Returns:
        list: Parsed intervals
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    all_intervals = []
    for line in tqdm(lines):
        line = line.strip()
        if line.startswith("Intervals ="):
            interval_part = line.split('=', 1)[1].strip()
            intervals = eval(interval_part)
            all_intervals.extend(intervals)
    return all_intervals

def process_fixed_intervals(fixed_intervals):  
    """
    Convert dataframe of fixed intervals to list of tuples.
    """
    processed_intervals = []
    for row in fixed_intervals.itertuples():
        processed_intervals.append((row.start, row.end))
    return processed_intervals

def optimize_intersection_check(list_1, list_2):
    """
    Count intersections between two interval lists.
    """
    intersections_count = 0 
    for coord1 in list_1:
        for coord2 in list_2:
            if is_intersect(coord1, coord2):
                intersections_count += 1
    return intersections_count

def process_file(se_files, names):
    """
    Process all super-enhancer files and calculate p-values.
    """
    data = {
        'HeLa':[], 'HCT116':[], 'MDA-MB-231':[],
        'MCF7':[], 'BT-549':[], 'LNCap':[]
    }
    
    # Process each chromosome
    for i, file in enumerate(z_dna_generations):
        combined_intervals = parse_intervals(file)
        chrom = 'chr' + str(i + 1)
        real_data = process_fixed_intervals(read_coordinates('ZDNABERT_hg38_generations.bed', chrom))
        
        # Calculate p-values for each cell line
        for se_file, name in zip(se_files, names):
            pval = p_value(real_data, combined_intervals, se_file, chrom)
            data[name].append(pval)
            print(data)
    return data

# File lists for different analyses
z_dna_generations = [f'whole_genome_generated_z_dna_on_chr{i}.txt' for i in range(1, 23)]
se_element_files = [
    'HeLa_SE_ele_hg38.bed', 'HCT116_SE_ele_hg38.bed', 
    'MDA-MB-231_SE_ele_hg38.bed', 'MCF7_SE_ele_hg38.bed',
    'BT_549_SE_ele_hg38.bed', 'LNCap_SE_ele_hg38.bed'
]
# te_files = [
#     'HeLa_SE_TE_hg38.bed', 'HCT116_SE_TE_hg38.bed', 
#     'MDA-MB-231_SE_TE_hg38.bed', 'MCF7_SE_TE_hg38.bed',
#     'BT_549_SE_TE_hg38.bed', 'LNCap_SE_TE_hg38.bed'
# ]
names = ['HeLa', 'HCT116', 'MDA-MB-231', 'MCF7', 'BT-549', 'LNCap']
    
# Execute analyses
data = process_file(se_element_files, names)
df = (pd.DataFrame(data)).T
df['result_fisher_method'] = df.apply(fisher_method, axis=1)
df.to_excel('results_fisher_method_se_elements.xlsx')

# data = process_file(te_files, names)
# df = (pd.DataFrame(data)).T
# df['result_fisher_method'] = df.apply(fisher_method, axis=1)
# df.to_excel('results_fisher_method_te.xlsx')


# ======================
# 3. VCF ANNOTATION
# ======================

"""
Bash script for annotating VCF files with BED regions.
Run with: bash annotate_vcfs.sh
"""

#!/bin/bash

# Configuration
INPUT_DIR="1000genomes_data"
ANNOTATION_BED="annotation.bed"
TMP_DIR="$_"
mkdir -p "$TMP_DIR"

process_vcf() {
    local vcf="$1"
    echo "Processing $vcf"
    
    # Setup temporary files
    local OUTPUT_FILE="${vcf%.gz}"
    local TABIX_FILE="$vcf.tbi"
    local BED_FILE="$TMP_DIR/$(basename "$vcf").$$.bed"
    local INTERSECT_FILE="$TMP_DIR/$(basename "$vcf").$$.intersect"
    local HEADER_FILE="$TMP_DIR/$(basename "$vcf").$$.header"
    local BODY_FILE="$TMP_DIR/$(basename "$vcf").$$.body"

    # Create index if missing
    [ ! -f "$TABIX_FILE" ] && tabix -p vcf "$vcf"

    # Extract header
    zcat "$vcf" | grep "^#" > "$HEADER_FILE"

    # Process by chromosome
    local CHROMS=$(cut -f1 "$ANNOTATION_BED" | sort -u)
    > "$BED_FILE"
    
    for chrom in $CHROMS; do
        local REGIONS=$(grep -w "^$chrom" "$ANNOTATION_BED" | awk '{print $1":"$2"-"$3}' | tr '\n' ' ')
        
        [ -n "$REGIONS" ] && \
        tabix "$vcf" $REGIONS | \
        awk '!/^#/ && $5 != "<NON_REF>" {print $1 "\t" $2 "\t" $2 "\t" $0}' >> "$BED_FILE"
    done

    # Intersect with annotation
    bedtools intersect -a "$BED_FILE" -b "$ANNOTATION_BED" -wa -wb > "$INTERSECT_FILE"

    # Reformat output
    awk 'BEGIN {FS=OFS="\t"} {
        printf "%s", $4;
        for (i=5; i<=NF; i++) printf "%s%s", OFS, $i;
        print "";
    }' "$INTERSECT_FILE" > "$BODY_FILE"

    # Combine and clean up
    cat "$HEADER_FILE" "$BODY_FILE" > "$OUTPUT_FILE"
    rm "$BED_FILE" "$INTERSECT_FILE" "$HEADER_FILE" "$BODY_FILE"
}

export -f process_vcf
export ANNOTATION_BED TMP_DIR

# Parallel processing
find "$INPUT_DIR" -name "*.vcf.gz" -print0 | \
  parallel -0 -j$(nproc) --progress --joblog "$TMP_DIR/parallel_joblog.txt" \
  "process_vcf {}"


# ======================
# 4. MUTATION CONCENTRATION
# ======================

def clean_dataframe(df):
    """
    Clean and validate mutation dataframe.
    
    Parameters:
        df: Input dataframe
        
    Returns:
        Cleaned dataframe with:
        - Stripped strings
        - Proper numeric types
        - Valid intervals
        - Removed duplicates
    """
    # Clean string columns
    str_cols = df.select_dtypes(include=['object']).columns
    for col in str_cols:
        df[col] = df[col].str.strip()
        df[col] = df[col].str.replace(r'\s+', ' ', regex=True)
    
    # Convert numeric columns
    df['INT_START'] = pd.to_numeric(df['INT_START'], errors='coerce')
    df['INT_END'] = pd.to_numeric(df['INT_END'], errors='coerce')
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    
    # Remove duplicates
    initial_rows = len(df)
    df = df.drop_duplicates()
    removed_rows = initial_rows - len(df)
    if removed_rows > 0:
        print(f"Removed {removed_rows} duplicate rows")
    
    # Validate intervals
    df = df[(df['INT_END'] > df['INT_START']) & 
            (df['INT_START'].notna()) & 
            (df['INT_END'].notna())]
    
    return df.reset_index(drop=True)

def calculate_mutation_concentration(df):
    """
    Calculate mutation concentration per genomic region.
    
    Parameters:
        df: Cleaned mutation dataframe
        
    Returns:
        Dataframe with:
        - SNP counts per region
        - Mutation concentrations
    """
    # Check for duplicates
    dup_check = df.duplicated(subset=['ID', 'CHROM', 'POS', 'GROUP'], keep=False)
    if dup_check.any():
        print("Warning: Duplicate entries found")
        display(df[dup_check].sort_values(['ID', 'CHROM', 'POS']).head())
    
    # Calculate SNP counts and concentrations
    snp_counts = (df.groupby(['ID', 'GROUP', 'INT_START', 'INT_END'])
                  ['POS'].nunique()
                  .reset_index(name='SNP_COUNT'))
    
    snp_counts['NEW_CONC'] = (snp_counts['SNP_COUNT'] / 
                             (snp_counts['INT_END'] - snp_counts['INT_START']))
    
    # Handle infinite values
    snp_counts['NEW_CONC'].replace([np.inf, -np.inf], np.nan, inplace=True)
    
    return snp_counts

# ==============================
# 5. SEQUENCE GENERATION FROM MUTATIONS
# ==============================

def generate_mutant_sequences(mutation_file, patient_data_file, reference_genome_file, 
                             homer_motifs_file, meme_motifs_file, output_folder="output_100_flank"):
    """
    Generate reference and mutant sequences with flanking regions from mutation data.
    
    Parameters:
        mutation_file (str): Path to mutation data Excel file
        patient_data_file (str): Path to patient data TSV file
        reference_genome_file (str): Path to reference genome FASTA
        homer_motifs_file (str): Path to HOMER motifs Excel file
        meme_motifs_file (str): Path to MEME motifs Excel file
        output_folder (str): Output directory for sequence files
        
    Returns:
        None (writes FASTA files to output directory)
    """
    from Bio import SeqIO
    
    # Mapping table to convert RefSeq identifiers to standard chromosome names
    refseq_to_chr = {
        "NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3",
        "NC_000004.12": "chr4", "NC_000005.10": "chr5", "NC_000006.12": "chr6",
        "NC_000007.14": "chr7", "NC_000008.11": "chr8", "NC_000009.12": "chr9",
        "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
        "NC_000013.11": "chr13", "NC_000014.9": "chr14", "NC_000015.10": "chr15",
        "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
        "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9": "chr21",
        "NC_000022.11": "chr22", "NC_000023.11": "chrX", "NC_000024.10": "chrY",
        "NC_012920.1": "chrMT"
    }

    # Load data files
    mutation_data = pd.read_excel(mutation_file)
    patient_data = pd.read_csv(patient_data_file, sep="\t")
    reference_genome = SeqIO.to_dict(SeqIO.parse(reference_genome_file, "fasta"))
    homer_motifs = pd.read_excel(homer_motifs_file)
    meme_motifs = pd.read_excel(meme_motifs_file)

    # Filter and normalize reference genome
    filtered_reference = {refseq_id: record for refseq_id, record in reference_genome.items() 
                         if refseq_id in refseq_to_chr}
    normalized_reference = {}
    for refseq_id, record in filtered_reference.items():
        if refseq_id in refseq_to_chr:
            normalized_reference[refseq_to_chr[refseq_id]] = record
        else:
            print(f"Warning: RefSeq ID {refseq_id} not found in mapping table. Skipping.")

    def extract_reference_sequence(chrom, start, end):
        """Extract reference sequence from genome"""
        if chrom in normalized_reference:
            return str(normalized_reference[chrom].seq[start-1:end])  # 1-based indexing
        raise ValueError(f"Chromosome {chrom} not found in reference genome.")

    def extract_sequence_with_flanks(chrom, start, end, flank_size=100):
        """Extract sequence with flanking regions"""
        if chrom in normalized_reference:
            seq_len = len(normalized_reference[chrom].seq)
            flanked_start = max(1, start - flank_size)
            flanked_end = min(seq_len, end + flank_size)
            return str(normalized_reference[chrom].seq[flanked_start-1:flanked_end])
        raise ValueError(f"Chromosome {chrom} not found in reference genome.")

    # Create output folder
    os.makedirs(output_folder, exist_ok=True)

    # Process each mutation
    for _, mutation in mutation_data.iterrows():
        mutation_type = mutation["type"]
        start = mutation["start"]
        end = mutation["end"]
        chrom = mutation["CHROM"]
        pos = mutation["pos"]

        # Create mutation-specific folder
        mutation_folder = os.path.join(output_folder, f"mutation_{mutation_type}_{chrom}_{start}_{end}_{pos}")
        os.makedirs(mutation_folder, exist_ok=True)

        try:
            # Extract reference sequences
            ref_seq = extract_reference_sequence(chrom, start, end)
            ref_seq_flanked = extract_sequence_with_flanks(chrom, start, end)
            
            # Write reference sequences
            with open(os.path.join(mutation_folder, "reference.fasta"), "w") as f:
                f.write(f">reference_{chrom}_{start}_{end}\n{ref_seq}\n")
            with open(os.path.join(mutation_folder, "reference_with_flanks.fasta"), "w") as f:
                f.write(f">reference_with_flanks_{chrom}_{start}_{end}\n{ref_seq_flanked}\n")

            # Find matching patient mutations
            matches = patient_data[(patient_data["CHROM"] == chrom) & (patient_data["POS"] == pos)]
            
            # Generate and write alternative sequences
            with open(os.path.join(mutation_folder, "merged_sequences.fasta"), "w") as merged:
                merged.write(f">reference_{chrom}_{start}_{end}\n{ref_seq}\n")
                merged.write(f">reference_with_flanks_{chrom}_{start}_{end}\n{ref_seq_flanked}\n")

                for _, row in matches.iterrows():
                    patient_id = row["ID"]
                    alt = row["ALT"]
                    all_ref = row["ALL"]
                    
                    # Generate mutant sequence
                    if mutation_type in ["snp", "sins", "sdel", "lins", "ldel"]:
                        alt_seq = ref_seq[:pos-start] + alt + ref_seq[pos-start+len(all_ref):]
                        alt_seq_flanked = extract_sequence_with_flanks(chrom, start, end)
                        alt_seq_flanked = (alt_seq_flanked[:100] + alt_seq + alt_seq_flanked[100+len(ref_seq):])
                        
                        # Write individual files
                        with open(os.path.join(mutation_folder, f"alternative_{patient_id}.fasta"), "w") as f:
                            f.write(f">{patient_id}_alternative_{chrom}_{start}_{end}\n{alt_seq}\n")
                        with open(os.path.join(mutation_folder, f"alternative_with_flanks_{patient_id}.fasta"), "w") as f:
                            f.write(f">{patient_id}_alternative_with_flanks_{chrom}_{start}_{end}\n{alt_seq_flanked}\n")
                        
                        # Add to merged file
                        merged.write(f">{patient_id}_alternative_with_flanks_{chrom}_{start}_{end}\n{alt_seq_flanked}\n")

                # Add motifs if present
                if mutation_type in homer_motifs.columns:
                    for motif in homer_motifs[mutation_type].dropna():
                        merged.write(f">homer_motif_{mutation_type}\n{motif}\n")
                if mutation_type in meme_motifs.columns:
                    for motif in meme_motifs[mutation_type].dropna():
                        merged.write(f">meme_motif_{mutation_type}\n{motif}\n")

        except ValueError as e:
            print(f"Error processing mutation {chrom}:{start}-{end}: {e}")
            continue
