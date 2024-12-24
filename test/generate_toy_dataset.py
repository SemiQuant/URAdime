import pysam
from Bio.Seq import Seq
import random
import pandas as pd
import os
from pathlib import Path

# Set fixed random seed for reproducibility
random.seed(42)

def create_test_primers():
    """Create test primers file with known sequences"""
    primers = {
        'Name': ['Primer1', 'Primer2', 'Primer3'],
        'Forward': [
            'TAATAAGCCCCCGTCACTGTTGGTTGT',
            'CCCAGGACGGGTTGGCCAGATGTG',
            'GCTTAGTGGCTCTTGGGCCGCGGTGCGTT'
        ],
        'Reverse': [
            'TTGTCCTTTTATCCGCTCACTT',
            'AAGCTTAAATGGGAAATACGCGGCCATAAG',
            'GAGAGCCAGCTGCGTTCGCTAATGTGAG'
        ],
        'Size': [200, 250, 400]  # Expected amplicon sizes
    }
    
    # Ensure test directory exists
    Path('test/data').mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(primers)
    df.to_csv('test/data/test_primers.tsv', sep='\t', index=False)
    return primers

def get_fixed_sequence(length, seed_multiplier=1):
    """Generate a fixed sequence of specified length"""
    # Use a combination of repeating patterns to ensure uniqueness
    random.seed(42 * seed_multiplier)  # Different multiplier for different calls
    return ''.join(random.choice('ACGT') for _ in range(length))

def generate_perfect_read(forward_primer, reverse_primer, size, read_num):
    """Generate a perfect read with correct orientation and size"""
    middle_length = size - len(forward_primer) - len(reverse_primer)
    middle = get_fixed_sequence(middle_length, read_num)
    sequence = forward_primer + middle + str(Seq(reverse_primer).reverse_complement())
    return sequence

def generate_wrong_size_read(forward_primer, reverse_primer, target_size, read_num):
    """Generate a read with correct primers but wrong size"""
    wrong_size = target_size + int(target_size * 0.2)  # 20% larger
    middle_length = wrong_size - len(forward_primer) - len(reverse_primer)
    middle = get_fixed_sequence(middle_length, read_num)
    sequence = forward_primer + middle + str(Seq(reverse_primer).reverse_complement())
    return sequence

def generate_single_terminal_match(primer, size, at_start=True, terminal_length=15, read_num=0):
    """Generate a read with a single terminal match"""
    # Calculate middle length to achieve target size
    middle_length = size - terminal_length
    fixed_seq = get_fixed_sequence(middle_length, read_num)
    
    # Only use part of the primer to ensure it's not detected as a full match
    partial = primer[:terminal_length]
    
    # Add multiple mismatches to prevent full match detection
    partial = partial[:-3] + 'NNN'
    
    if at_start:
        # For start matches, use forward primer sequence
        return partial + fixed_seq
    else:
        # For end matches, use reverse complement
        return fixed_seq + str(Seq(partial).reverse_complement())

def generate_paired_terminal_match(forward_primer, reverse_primer, size, terminal_length=15, read_num=0):
    """Generate a read with terminal matches at both ends"""
    # Calculate middle length to achieve target size
    middle_length = size - (2 * terminal_length)
    middle = get_fixed_sequence(middle_length, read_num)
    
    # Only use part of the primers to ensure they're not detected as full matches
    partial_fwd = forward_primer[:terminal_length]
    partial_rev = reverse_primer[:terminal_length]
    
    # Add multiple mismatches to prevent full match detection
    partial_fwd = partial_fwd[:-3] + 'NNN'
    partial_rev = partial_rev[:-3] + 'NNN'
    
    # Ensure correct orientation: forward at start, reverse complement at end
    return partial_fwd + middle + str(Seq(partial_rev).reverse_complement())

def generate_hybrid_match(forward_primer, reverse_primer, size, use_forward=True, terminal_length=15, read_num=0):
    """Generate a read with one full primer and one terminal match"""
    if use_forward:
        # Full forward primer + terminal match of reverse primer
        middle_length = size - len(forward_primer) - terminal_length
        middle = get_fixed_sequence(middle_length, read_num)
        partial_rev = reverse_primer[:terminal_length]
        # Add multiple mismatches to terminal match
        partial_rev = partial_rev[:-3] + 'NNN'
        sequence = forward_primer + middle + str(Seq(partial_rev).reverse_complement())
    else:
        # Terminal match of forward primer + full reverse primer
        middle_length = size - terminal_length - len(reverse_primer)
        middle = get_fixed_sequence(middle_length, read_num)
        partial_fwd = forward_primer[:terminal_length]
        # Add multiple mismatches to terminal match
        partial_fwd = partial_fwd[:-3] + 'NNN'
        sequence = partial_fwd + middle + str(Seq(reverse_primer).reverse_complement())
    
    return sequence

def create_test_bam():
    """Create test BAM file with a balanced distribution of reads"""
    header = {'HD': {'VN': '1.0'},
             'SQ': [{'LN': 1000, 'SN': 'test_ref'}]}
    
    primers = create_test_primers()
    read_count = 0
    
    # Ensure test directory exists
    Path('test/data').mkdir(parents=True, exist_ok=True)
    
    with pysam.AlignmentFile('test/data/test_reads.bam', 'wb', header=header) as outf:
        # Generate different types of reads for each primer pair
        for i, (name, fwd, rev, size) in enumerate(zip(primers['Name'], 
                                                      primers['Forward'], 
                                                      primers['Reverse'], 
                                                      primers['Size'])):
        
            # 1. Perfect matches (40 reads per primer)
            for j in range(40):
                sequence = generate_perfect_read(fwd, rev, size, read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"perfect_match_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 2. Wrong size matches (40 reads per primer)
            for j in range(40):
                sequence = generate_wrong_size_read(fwd, rev, size, read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"wrong_size_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 3. Single terminal matches (20 correct size, 20 wrong size per primer)
            for j in range(20):
                # Start terminal match with correct size
                sequence = generate_single_terminal_match(fwd, size, at_start=True, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_start_correct_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

                # End terminal match with correct size
                sequence = generate_single_terminal_match(rev, size, at_start=False, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_end_correct_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1
                
                # Start terminal match with wrong size
                wrong_size = int(size * 1.2)  # 20% larger
                sequence = generate_single_terminal_match(fwd, wrong_size, at_start=True, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_start_wrong_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

                # End terminal match with wrong size
                sequence = generate_single_terminal_match(rev, wrong_size, at_start=False, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_end_wrong_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 4. Paired terminal matches (20 correct size, 20 wrong size per primer)
            for j in range(20):
                # Correct size
                sequence = generate_paired_terminal_match(fwd, rev, size, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"paired_terminal_correct_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1
                
                # Wrong size
                wrong_size = int(size * 1.2)  # 20% larger
                sequence = generate_paired_terminal_match(fwd, rev, wrong_size, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"paired_terminal_wrong_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 5. Single-end primers (20 forward, 20 reverse per primer)
            for j in range(20):
                # Forward primer only
                fixed_seq = get_fixed_sequence(100, read_count)
                sequence = fwd + fixed_seq
                a = pysam.AlignedSegment()
                a.query_name = f"single_end_forward_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

                # Reverse primer only
                fixed_seq = get_fixed_sequence(100, read_count)
                sequence = fixed_seq + str(Seq(rev).reverse_complement())
                a = pysam.AlignedSegment()
                a.query_name = f"single_end_reverse_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 6. Hybrid matches (20 correct size, 20 wrong size per primer)
            for j in range(20):
                # Correct size with full forward + terminal reverse
                sequence = generate_hybrid_match(fwd, rev, size, use_forward=True, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"hybrid_forward_full_correct_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1
                
                # Wrong size with terminal forward + full reverse
                wrong_size = int(size * 1.2)  # 20% larger
                sequence = generate_hybrid_match(fwd, rev, wrong_size, use_forward=False, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"hybrid_reverse_full_wrong_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 7. Mismatched primer pairs (40 reads per primer)
            for j in range(40):
                # Use primers from different pairs
                other_fwd = primers['Forward'][(i + 1) % 3]
                other_rev = primers['Reverse'][(i + 2) % 3]
                middle = get_fixed_sequence(100, read_count)
                sequence = other_fwd + middle + str(Seq(other_rev).reverse_complement())
                a = pysam.AlignedSegment()
                a.query_name = f"mismatched_pair_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 8. No primers (40 reads per primer set)
            for j in range(40):
                sequence = get_fixed_sequence(100, read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"no_primers_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 9. Multi-primer pairs (40 reads per primer)
            for j in range(40):
                # Create sequence with multiple primers at one end
                other_fwd = primers['Forward'][(i + 1) % 3]
                middle = get_fixed_sequence(50, read_count)
                sequence = fwd + other_fwd + middle + str(Seq(rev).reverse_complement())
                a = pysam.AlignedSegment()
                a.query_name = f"multi_primer_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

    # Sort and index BAM file
    pysam.sort('-o', 'test/data/test_data.bam', 'test/data/test_reads.bam')
    os.remove('test/data/test_reads.bam')
    pysam.index('test/data/test_data.bam')
    
    print(f"Total reads created: {read_count}")

if __name__ == "__main__":
    create_test_bam()