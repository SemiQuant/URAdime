# URAdime

URAdime (Universal Read Analysis of DIMErs) is a Python package for analyzing primer sequences in sequencing data to identify dimers and chimeras.

## Installation

```bash
pip install uradime
```

## Usage

URAdime can be used both as a command-line tool and as a Python package.

### Command Line Interface

```bash
# Basic usage
uradime -b input.bam -p primers.tsv -o results/my_analysis

# Full options
uradime \
    -b input.bam \                    # Input BAM file
    -p primers.tsv \                  # Primer file (tab-separated)
    -o results/my_analysis \          # Output prefix
    -t 8 \                            # Number of threads
    -m 1000 \                         # Maximum reads to process (0 for all)
    -c 100 \                          # Chunk size for parallel processing
    -u \                              # Process only unaligned reads
    --max-distance 2 \                # Maximum Levenshtein distance for matching
    --unaligned-only \                # only check the unaligned reads  
    --window-size 20 \                # Allowed padding on the 5' ends of the reads, sometime needs to be very big due to universal tails etc. setting this parameter too large can cause unexpected results
    --ignore-amplicon-size \          # Usefull if short read sequecing like Illumina where the paired read length is not the size of the actual amplicon
    --check-termini \                 # Turn off check for partial matches at read termini
    --terminus-length 14 \            # Length of terminus to check for partial matches
    --overlap-threshold 0.8 \         # Minimum fraction of overlap required to consider primers as overlapping (0.0-1.0), this is added for hissPCR support
    --downsample 5.0 \                # Percentage of reads to randomly sample from the BAM file (0.1-100.0)
    -v                                # Verbose output
```



### Python Package

```python
from uradime import bam_to_fasta_parallel, create_analysis_summary, load_primers

# Load and analyze BAM file
result_df = bam_to_fasta_parallel(
    bam_path="your_file.bam",
    primer_file="primers.tsv",
    num_threads=4
)

# Load primers for analysis
primers_df, _ = load_primers("primers.tsv")

# Create analysis summary
summary_df, matched_pairs, mismatched_pairs = create_analysis_summary(result_df, primers_df)
```

## Input Files

### Primer File Format (TSV)
The primer file should be tab-separated with the following columns:
- Name: Primer pair name
- Forward: Forward primer sequence
- Reverse: Reverse primer sequence
- Size: Expected amplicon size

Example:
```
Name    Forward             Reverse             Size
Pair1   ATCGATCGATCG       TAGCTAGCTAGC       100
Pair2   GCTAGCTAGCTA       CGATTCGATCGA       150
```

## Output Files

The tool generates several CSV files with the analysis results:
- `*_summary.csv`: Overall analysis summary
- `*_matched_pairs.csv`: Reads with matching primer pairs
- `*_mismatched_pairs.csv`: Reads with mismatched primer pairs
- `*_wrong_size_pairs.csv`: Reads with correct primer pairs but wrong size


## Requirements

- Python ≥3.7
- pysam
- pandas
- biopython
- python-Levenshtein
- tqdm
- numpy

## License

This project is licensed under GNU GPL.