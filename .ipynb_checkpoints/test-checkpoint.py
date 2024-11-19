import pysam
import pandas as pd
from Bio.Seq import Seq
import Levenshtein as lev

def load_primers(primer_file):
    """Load and prepare primers dataframe"""
    primers_df = pd.read_csv(primer_file, sep="\t")
    primers_df = primers_df.dropna(subset=['Forward', 'Reverse'])
    longest_primer_length = max(
        primers_df['Forward'].apply(len).max(), 
        primers_df['Reverse'].apply(len).max()
    )
    return primers_df, longest_primer_length
def is_match(seq1, seq2, max_distance=2):
    """
    Check for approximate match using Levenshtein distance.
    seq1: longer sequence to search in
    seq2: primer sequence to find
    max_distance: maximum allowed edit distance
    """
    # Handle potential None or empty sequences
    if not seq1 or not seq2:
        return False
    
    try:
        # Slide the primer (seq2) across the sequence (seq1)
        for i in range(len(seq1) - len(seq2) + 1):
            window = seq1[i:i+len(seq2)]
            if len(window) == len(seq2):  # Ensure we have a full window
                distance = lev.distance(str(window), str(seq2))
                if distance <= max_distance:
                    return True
    except:
        return False
    return False

def find_primers_in_region(sequence, primers_df, window_size=100, max_distance=2):
    """Find primers in a given sequence region"""
    primers_found = []
    
    for _, primer in primers_df.iterrows():
        forward_primer = primer['Forward']
        reverse_primer = primer['Reverse']
        reverse_complement_forward = str(Seq(forward_primer).reverse_complement())
        reverse_complement_reverse = str(Seq(reverse_primer).reverse_complement())
        
        # Check each primer against the entire sequence
        if is_match(sequence, forward_primer, max_distance):
            primers_found.append(f"{primer['Name']}_Forward")
            
        if is_match(sequence, reverse_primer, max_distance):
            primers_found.append(f"{primer['Name']}_Reverse")
            
        if is_match(sequence, reverse_complement_forward, max_distance):
            primers_found.append(f"{primer['Name']}_ForwardComp")
            
        if is_match(sequence, reverse_complement_reverse, max_distance):
            primers_found.append(f"{primer['Name']}_ReverseComp")
    
    return list(set(primers_found))  # Remove duplicates
    
def bam_to_fasta(bam_path, primer_file, unaligned_only=False, max_reads=200):
    """Process BAM file and find primers in reads"""
    # Load primers
    primers_df, longest_primer_length = load_primers(primer_file)
    
    # Open BAM file
    try:
        bam_file = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"Error opening BAM file: {e}")
        return pd.DataFrame()

    data = []
    search_window = 100 + longest_primer_length  # Increased search window size
    
    reads_processed = 0
    
    for read in bam_file.fetch(until_eof=True):
        if unaligned_only and not read.is_unmapped:
            continue
        if read.query_sequence is None:
            continue

        reads_processed += 1
        
        if reads_processed > max_reads:
            break
        
        read_sequence = read.query_sequence
        read_length = len(read_sequence)
        
        # print(f"\nProcessing read {reads_processed}:")
        # print(f"Read name: {read.query_name}")
        # print(f"Read length: {read_length}")
        
        # Define search regions with bounds checking
        start_region = read_sequence[:min(search_window, read_length)]
        end_region = read_sequence[max(0, read_length - search_window):]
        
        # Find primers in both regions
        start_primers_found = find_primers_in_region(start_region, primers_df, window_size=search_window, max_distance=2)
        end_primers_found = find_primers_in_region(end_region, primers_df, window_size=search_window, max_distance=2)
        
        data.append({
            'Read_Name': read.query_name,
            'Start_Primers': ', '.join(start_primers_found) if start_primers_found else 'None',
            'End_Primers': ', '.join(end_primers_found) if end_primers_found else 'None',
            'Read_Length': read_length
        })
    
    bam_file.close()
    
    result_df = pd.DataFrame(data)
    return result_df

def main():
    try:
        result_df = bam_to_fasta(
            bam_path="test.bam",
            primer_file="primers.tsv",
            unaligned_only=False
        )
        print("\nResults summary:")
        print(f"Total reads processed: {len(result_df)}")
        print("\nDetailed results:")
        print(result_df)
        
        # Print statistics about primer findings
        print("\nPrimer statistics:")
        start_primers_found = result_df['Start_Primers'].str.count(',').sum() + len(result_df[result_df['Start_Primers'] != 'None'])
        end_primers_found = result_df['End_Primers'].str.count(',').sum() + len(result_df[result_df['End_Primers'] != 'None'])
        print(f"Total primers found at start: {start_primers_found}")
        print(f"Total primers found at end: {end_primers_found}")
        
    except Exception as e:
        print(f"Error in main execution: {e}")

if __name__ == "__main__":
    main()


def split_primer_results(result_df):
    """
    Split results into three dataframes:
    1. single_primer_df: Exactly one primer at each end
    2. multiple_primer_df: More than one primer at either end
    3. no_primer_df: No primers at either end
    """
    
    def count_primers(row):
        """Helper function to count primers in Start_Primers and End_Primers"""
        start_count = 0 if row['Start_Primers'] == 'None' else row['Start_Primers'].count(',') + 1
        end_count = 0 if row['End_Primers'] == 'None' else row['End_Primers'].count(',') + 1
        return start_count, end_count

    # Create mask for each category
    single_primer_mask = result_df.apply(
        lambda row: count_primers(row) == (1, 1), 
        axis=1
    )
    
    multiple_primer_mask = result_df.apply(
        lambda row: any(count > 1 for count in count_primers(row)), 
        axis=1
    )
    
    no_primer_mask = result_df.apply(
        lambda row: count_primers(row) == (0, 0), 
        axis=1
    )
    
    # Split into three dataframes
    single_primer_df = result_df[single_primer_mask].copy()
    multiple_primer_df = result_df[multiple_primer_mask].copy()
    no_primer_df = result_df[no_primer_mask].copy()
    
    # Add summary information
    print("\nSummary of split results:")
    print(f"Reads with single primer at each end: {len(single_primer_df)}")
    print(f"Reads with multiple primers at either end: {len(multiple_primer_df)}")
    print(f"Reads with no primers: {len(no_primer_df)}")
    
    print("\nSingle primer reads:")
    print(single_primer_df)
    
    print("\nMultiple primer reads:")
    print(multiple_primer_df)
    
    print("\nNo primer reads:")
    print(no_primer_df)
    
    return single_primer_df, multiple_primer_df, no_primer_df

# Modify the main function to include the splitting
def main():
    try:
        result_df = bam_to_fasta(
            bam_path="test_sub.bam",
            primer_file="primers.tsv",
            unaligned_only=False
        )
        print("\nInitial results summary:")
        print(f"Total reads processed: {len(result_df)}")
        print("\nAll results:")
        print(result_df)
        
        # Split results into categories
        single_primer_df, multiple_primer_df, no_primer_df = split_primer_results(result_df)
        
        # Save results to CSV files (optional)
        single_primer_df.to_csv('single_primer_reads.csv', index=False)
        multiple_primer_df.to_csv('multiple_primer_reads.csv', index=False)
        no_primer_df.to_csv('no_primer_reads.csv', index=False)
        
    except Exception as e:
        print(f"Error in main execution: {e}")

if __name__ == "__main__":
    main()

def analyze_primer_pairs(single_primer_df, primers_df):
    """
    Analyze single primer matches for expected F/R pairs and size compliance
    """
    # Helper function to extract primer name without orientation
    def get_base_primer_name(primer_str):
        if primer_str == 'None':
            return None
        return primer_str.split('_')[0]  # Remove the _Forward or _Reverse suffix
    
    # Helper function to get primer orientation
    def get_primer_orientation(primer_str):
        if primer_str == 'None':
            return None
        return primer_str.split('_')[1]  # Get Forward or Reverse part
    
    # Add columns for base primer names and orientations
    single_primer_df['Start_Primer_Name'] = single_primer_df['Start_Primers'].apply(get_base_primer_name)
    single_primer_df['End_Primer_Name'] = single_primer_df['End_Primers'].apply(get_base_primer_name)
    single_primer_df['Start_Orientation'] = single_primer_df['Start_Primers'].apply(get_primer_orientation)
    single_primer_df['End_Orientation'] = single_primer_df['End_Primers'].apply(get_primer_orientation)
    
    # Find matching pairs (same primer name at both ends)
    matching_pairs_df = single_primer_df[
        single_primer_df['Start_Primer_Name'] == single_primer_df['End_Primer_Name']
    ].copy()
    
    # Check correct orientation (one Forward, one Reverse)
    correct_orientation_df = matching_pairs_df[
        ((matching_pairs_df['Start_Orientation'].str.contains('Forward') & 
          matching_pairs_df['End_Orientation'].str.contains('Reverse')) |
         (matching_pairs_df['Start_Orientation'].str.contains('Reverse') & 
          matching_pairs_df['End_Orientation'].str.contains('Forward')))
    ].copy()
    
    # Add expected size information from primers_df
    primer_sizes = primers_df.set_index('Name')['Size'].to_dict()
    correct_orientation_df['Expected_Size'] = correct_orientation_df['Start_Primer_Name'].map(primer_sizes)
    
    # Calculate size compliance (within 10% of expected)
    def is_size_compliant(row):
        expected = row['Expected_Size']
        actual = row['Read_Length']
        if pd.isna(expected):
            return False
        tolerance = expected * 0.10  # 10% tolerance
        return abs(actual - expected) <= tolerance
    
    correct_orientation_df['Size_Compliant'] = correct_orientation_df.apply(is_size_compliant, axis=1)
    
    # Print summary statistics
    print("\nPrimer Pair Analysis Summary:")
    print(f"Total single primer matches: {len(single_primer_df)}")
    print(f"Matching primer pairs (same primer at both ends): {len(matching_pairs_df)}")
    print(f"Correct orientation pairs: {len(correct_orientation_df)}")
    print(f"Size compliant pairs: {len(correct_orientation_df[correct_orientation_df['Size_Compliant']])}")
    
    print("\nDetailed results for correct orientation pairs:")
    print(correct_orientation_df[[
        'Read_Name', 'Start_Primers', 'End_Primers', 
        'Read_Length', 'Expected_Size', 'Size_Compliant'
    ]])
    
    # Save detailed results
    correct_orientation_df.to_csv('primer_pair_analysis.csv', index=False)
    
    return correct_orientation_df

# Modify the main function to include the primer pair analysis
def main():
    try:
        result_df = bam_to_fasta(
            bam_path="test_sub.bam",
            primer_file="primers.tsv",
            unaligned_only=False
        )
        
        # Load primers for size information
        primers_df, _ = load_primers("primers.tsv")
        
        # Split results into categories
        single_primer_df, multiple_primer_df, no_primer_df = split_primer_results(result_df)
        
        # Analyze primer pairs in single_primer_df
        correct_pairs_df = analyze_primer_pairs(single_primer_df, primers_df)
        
        # Additional statistics about non-matching pairs
        print("\nBreakdown of non-matching pairs:")
        non_matching = single_primer_df[
            single_primer_df['Start_Primer_Name'] != single_primer_df['End_Primer_Name']
        ]
        print("\nReads with different primers at each end:")
        print(non_matching[['Read_Name', 'Start_Primers', 'End_Primers', 'Read_Length']])
        
    except Exception as e:
        print(f"Error in main execution: {e}")

if __name__ == "__main__":
    main()