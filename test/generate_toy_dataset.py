import pysam
from Bio.Seq import Seq
import random
import pandas as pd
import os

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
    
    df = pd.DataFrame(primers)
    df.to_csv('test_primers.tsv', sep='\t', index=False)
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

def generate_single_terminal_match(primer, at_start=True, terminal_length=15, read_num=0):
    """Generate a read with a single terminal match"""
    buffer = "GGGGGG"
    fixed_seq = get_fixed_sequence(100, read_num)
    partial = primer[:terminal_length]
    if at_start:
        return partial + buffer + fixed_seq
    else:
        # For end matches, use reverse complement if it's meant to be at the end
        return fixed_seq + buffer + str(Seq(partial).reverse_complement())

def generate_paired_terminal_match(forward_primer, reverse_primer, terminal_length=15, read_num=0):
    """Generate a read with terminal matches at both ends"""
    buffer = "GGGGGG"
    middle = get_fixed_sequence(100, read_num)
    partial_fwd = forward_primer[:terminal_length]
    partial_rev = reverse_primer[:terminal_length]
    return partial_fwd + buffer + middle + buffer + str(Seq(partial_rev).reverse_complement())

def generate_hybrid_match(forward_primer, reverse_primer, use_forward=True, terminal_length=15, read_num=0):
    """Generate a read with one full primer and one terminal match"""
    buffer = "GGGGGG"
    middle = get_fixed_sequence(100, read_num)
    
    if use_forward:
        # Full forward primer + terminal match of reverse primer
        partial_rev = reverse_primer[:terminal_length]
        sequence = forward_primer + buffer + middle + buffer + str(Seq(partial_rev).reverse_complement())
    else:
        # Terminal match of forward primer + full reverse primer
        partial_fwd = forward_primer[:terminal_length]
        sequence = partial_fwd + buffer + middle + buffer + str(Seq(reverse_primer).reverse_complement())
    
    return sequence

def create_test_bam():
    """Create test BAM file with a balanced distribution of reads"""
    header = {'HD': {'VN': '1.0'},
             'SQ': [{'LN': 1000, 'SN': 'test_ref'}]}
    
    primers = create_test_primers()
    read_count = 0
    buffer = "GGGGGG"
    
    with pysam.AlignmentFile('test_reads.bam', 'wb', header=header) as outf:
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

            # 3. Single terminal matches (20 at start, 20 at end per primer)
            for j in range(20):
                # Start terminal match
                sequence = generate_single_terminal_match(fwd, at_start=True, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_start_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

                # End terminal match
                sequence = generate_single_terminal_match(rev, at_start=False, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"single_terminal_end_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 4. Paired terminal matches (40 reads per primer)
            for j in range(40):
                sequence = generate_paired_terminal_match(fwd, rev, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"paired_terminal_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 5. Single-end primers (20 forward, 20 reverse per primer)
            for j in range(20):
                # Forward primer only
                fixed_seq = get_fixed_sequence(100, read_count)
                sequence = fwd + buffer + fixed_seq
                a = pysam.AlignedSegment()
                a.query_name = f"single_end_forward_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

                # Reverse primer only
                fixed_seq = get_fixed_sequence(100, read_count)
                sequence = fixed_seq + buffer + str(Seq(rev).reverse_complement())
                a = pysam.AlignedSegment()
                a.query_name = f"single_end_reverse_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 6. Mismatched primer pairs (40 reads per primer)
            for j in range(40):
                # Use primers from different pairs
                other_fwd = primers['Forward'][(i + 1) % 3]
                other_rev = primers['Reverse'][(i + 2) % 3]
                middle = get_fixed_sequence(100, read_count)
                sequence = other_fwd + buffer + middle + buffer + str(Seq(other_rev).reverse_complement())
                a = pysam.AlignedSegment()
                a.query_name = f"mismatched_pair_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 7. No primers (40 reads per primer set)
            for j in range(40):
                sequence = get_fixed_sequence(100, read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"no_primers_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

            # 8. Hybrid matches (20 with full forward + terminal reverse, 20 with terminal forward + full reverse per primer)
            for j in range(20):
                # Full forward primer + terminal reverse
                sequence = generate_hybrid_match(fwd, rev, use_forward=True, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"hybrid_forward_full_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1
                
                # Terminal forward + full reverse
                sequence = generate_hybrid_match(fwd, rev, use_forward=False, read_num=read_count)
                a = pysam.AlignedSegment()
                a.query_name = f"hybrid_reverse_full_{name}_{j}"
                a.query_sequence = sequence
                a.flag = 4
                a.is_unmapped = True
                outf.write(a)
                read_count += 1

    # Sort and index BAM file
    pysam.sort('-o', 'test_data.bam', 'test_reads.bam')
    os.remove('test_reads.bam')
    pysam.index('test_data.bam')
    
    print(f"Total reads created: {read_count}")

if __name__ == "__main__":
    create_test_bam()