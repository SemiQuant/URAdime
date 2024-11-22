#!/usr/bin/env python3
"""
URAdime (Universal Read Analysis of DIMErs)
A tool for analyzing primer dimers and other artifacts in sequencing data.
This script processes BAM files to identify primer sequences at read ends
and analyzes their relationships and orientations.
"""

import argparse
import pysam
import pandas as pd
from Bio.Seq import Seq
import Levenshtein as lev
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional
import numpy as np
from tqdm import tqdm
import os
import sys
import random

def print_banner():
    """Print a stylized banner."""
    banner = """
══════════════════════════════════════════════════════════════════════════════
                           URAdime v{:^7}                
                Universal Read Analysis of DIMErs    
══════════════════════════════════════════════════════════════════════════════""".format("0.1.6")
    print(banner)

def load_primers(primer_file):
    """
    Load and prepare primers from a tab-separated file.
    Automatically replaces commas and semicolons in primer names with underscores.
    
    Args:
        primer_file (str): Path to tab-separated primer file containing Name, Forward, Reverse, and Size columns
        
    Returns:
        tuple: (primers_df, longest_primer_length)
            - primers_df: DataFrame containing primer information
            - longest_primer_length: Length of the longest primer sequence
            
    Raises:
        FileNotFoundError: If primer file doesn't exist
    """
    if not os.path.exists(primer_file):
        raise FileNotFoundError(f"Primer file not found: {primer_file}")
        
    primers_df = pd.read_csv(primer_file, sep="\t")
    primers_df = primers_df.dropna(subset=['Forward', 'Reverse'])
    
    # Clean primer names by replacing commas and semicolons with underscores
    primers_df['Name'] = primers_df['Name'].apply(lambda x: str(x).replace(',', '_').replace(';', '_'))
    
    # Check for duplicate names after cleaning
    duplicate_names = primers_df['Name'].duplicated()
    if duplicate_names.any():
        duplicate_list = primers_df[duplicate_names]['Name'].tolist()
        raise ValueError(f"Duplicate primer names found after cleaning: {', '.join(duplicate_list)}")
    
    longest_primer_length = max(
        primers_df['Forward'].apply(len).max(), 
        primers_df['Reverse'].apply(len).max()
    )
    return primers_df, longest_primer_length

def is_match(seq1, seq2, max_distance):
    """
    Check for approximate match between two sequences using Levenshtein distance.
    Handles case-insensitive comparison and treats N's as potential matches.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        max_distance (int): Maximum allowed Levenshtein distance
        
    Returns:
        bool: True if sequences match within specified distance
    """
    if not seq1 or not seq2:
        return False
    
    try:
        # Convert to uppercase and remove N's for comparing lengths
        seq1 = str(seq1).upper()
        seq2 = str(seq2).upper()
        
        for i in range(len(seq1) - len(seq2) + 1):
            window = seq1[i:i+len(seq2)]
            if len(window) == len(seq2):
                # Calculate distance considering N's as potential matches
                distance = 0
                for w, s in zip(window, seq2):
                    if w != s and w != 'N' and s != 'N':
                        distance += 1
                        if distance > max_distance:
                            break
                
                if distance <= max_distance:
                    return True
    except:
        return False
    return False

def check_terminal_match(sequence, primer, terminus_length=15, max_distance=2):
    """
    Check for partial matches at sequence termini, including reverse complements.
    
    Args:
        sequence (str): DNA sequence to search
        primer (str): Primer sequence to look for
        terminus_length (int): Minimum length of terminus to check
        max_distance (int): Maximum allowed Levenshtein distance
    
    Returns:
        tuple: (bool, int) - (whether match found, length of longest match found)
    """
    if len(sequence) < terminus_length or len(primer) < terminus_length:
        return False, 0
        
    primer_rc = str(Seq(primer).reverse_complement())
    best_match_length = 0
    found_match = False
    
    # Check all possible combinations:
    # 1. Start of sequence vs start of primer
    # 2. Start of sequence vs RC of end of primer
    # 3. End of sequence vs RC of start of primer
    # 4. End of sequence vs end of primer
    
    combinations = [
        (sequence[:terminus_length], primer[:terminus_length]),  # Start vs Start
        (sequence[:terminus_length], primer_rc[-terminus_length:]),  # Start vs RC End
        (sequence[-terminus_length:], primer_rc[:terminus_length]),  # End vs RC Start
        (sequence[-terminus_length:], primer[-terminus_length:])  # End vs End
    ]
    
    for seq_part, primer_part in combinations:
        # Try progressively larger windows
        for window_size in range(terminus_length, min(len(sequence), len(primer)) + 1):
            if len(seq_part) >= window_size and len(primer_part) >= window_size:
                if is_match(seq_part[:window_size], primer_part[:window_size], max_distance):
                    found_match = True
                    best_match_length = max(best_match_length, window_size)
                else:
                    break  # Stop increasing window size if no match found
    
    return found_match, best_match_length

def find_primers_in_region(sequence, primers_df, window_size=20, max_distance=2, check_termini=True, terminus_length=15, overlap_threshold=0.8):
    """
    Search for primer sequences within a given region of DNA sequence.
    Handles overlapping primers by selecting the best match based on binding and orientation.
    
    Args:
        sequence: DNA sequence to search
        primers_df: DataFrame containing primer information
        window_size: Size of window to search
        max_distance: Maximum allowed Levenshtein distance
        check_termini: Whether to check for partial matches at termini
        terminus_length: Length of terminus to check
        overlap_threshold: Minimum overlap fraction to consider primers as overlapping
        
    Returns:
        dict: Dictionary containing full matches and terminal matches
    """
    full_matches = []
    terminal_matches = []
    
    # Store detailed match information for filtering overlaps
    match_details = []
    
    # First pass: find all possible matches with their positions
    for _, primer in primers_df.iterrows():
        forward_primer = primer['Forward']
        reverse_primer = primer['Reverse']
        
        # Check all possible orientations
        primer_checks = [
            (forward_primer, 'Forward'),
            (reverse_primer, 'Reverse'),
            (str(Seq(forward_primer).reverse_complement()), 'ForwardComp'),
            (str(Seq(reverse_primer).reverse_complement()), 'ReverseComp')
        ]
        
        for primer_seq, orientation in primer_checks:
            # Search for matches within the window
            for i in range(min(window_size, len(sequence) - len(primer_seq) + 1)):
                window = sequence[i:i+len(primer_seq)]
                if len(window) == len(primer_seq):
                    # Calculate mismatches considering N's
                    mismatches = sum(1 for w, p in zip(window, primer_seq) 
                                   if w != p and w != 'N' and p != 'N')
                    
                    if mismatches <= max_distance:
                        match_details.append({
                            'name': f"{primer['Name']}_{orientation}",
                            'start': i,
                            'end': i + len(primer_seq),
                            'length': len(primer_seq),
                            'mismatches': mismatches,
                            'sequence': primer_seq,
                            'is_correct_orientation': orientation in ['Forward', 'ReverseComp']
                        })
    
    # Filter overlapping matches
    if match_details:
        # Sort by start position and then by length (longer primers first)
        match_details.sort(key=lambda x: (x['start'], -x['length']))
        
        filtered_matches = []
        i = 0
        while i < len(match_details):
            current_match = match_details[i]
            overlapping_group = [current_match]
            
            # Find all matches that overlap with current match
            j = i + 1
            while j < len(match_details):
                next_match = match_details[j]
                
                # Calculate overlap
                overlap_start = max(current_match['start'], next_match['start'])
                overlap_end = min(current_match['end'], next_match['end'])
                overlap_length = max(0, overlap_end - overlap_start)
                
                # Calculate overlap fraction relative to shorter match
                min_length = min(current_match['length'], next_match['length'])
                overlap_fraction = overlap_length / min_length
                
                if overlap_fraction >= overlap_threshold:
                    overlapping_group.append(next_match)
                    j += 1
                else:
                    break
                    
            # Choose the best match from overlapping group
            if len(overlapping_group) > 1:
                # Score each match in the group
                best_match = max(overlapping_group, key=lambda x: (
                    -x['mismatches'],  # Fewer mismatches is better
                    x['is_correct_orientation'],  # Correct orientation is better
                    x['length'],  # Longer length is better
                    -x['start']  # Earlier start position is better
                ))
                filtered_matches.append(best_match)
            else:
                filtered_matches.append(current_match)
            
            i = j
        
        # Add filtered matches to results
        full_matches = [match['name'] for match in filtered_matches]
    
    # Handle terminal matches if needed and no full matches were found
    if check_termini and not full_matches:
        for _, primer in primers_df.iterrows():
            # Check forward primer terminal matches
            fwd_match_found, fwd_match_length = check_terminal_match(
                sequence, primer['Forward'], terminus_length, max_distance
            )
            if fwd_match_found:
                terminal_matches.append(f"{primer['Name']}_Forward_Terminal_{fwd_match_length}bp")
            
            # Check reverse primer terminal matches
            rev_match_found, rev_match_length = check_terminal_match(
                sequence, primer['Reverse'], terminus_length, max_distance
            )
            if rev_match_found:
                terminal_matches.append(f"{primer['Name']}_Reverse_Terminal_{rev_match_length}bp")
    
    return {
        'full_matches': list(set(full_matches)),
        'terminal_matches': list(set(terminal_matches))
    }

def bam_to_fasta_parallel(bam_path: str, primer_file: str, window_size: int = 20, 
                         unaligned_only: bool = False, max_reads: int = 200, 
                         num_threads: int = 4, chunk_size: int = 50, 
                         downsample_percentage: float = 100.0,
                         max_distance: int = 2,
                         overlap_threshold: float = 0.8) -> pd.DataFrame:
    """
    Process BAM file and find primers in reads using multiple threads.
    
    Args:
        ... (existing args)
        overlap_threshold: Minimum fraction of overlap required to consider primers as overlapping (0.0-1.0)
    """
    # Load primers
    primers_df, _ = load_primers(primer_file)
    
    print(f"Loading BAM file: {bam_path}")
    
    # Perform downsampling
    all_reads = downsample_reads(bam_path, downsample_percentage, max_reads)
    
    if not all_reads:
        print("No reads selected after downsampling")
        return pd.DataFrame()
    
    print(f"Processing {len(all_reads)} reads with {num_threads} threads...")
    
    chunks = [all_reads[i:i + chunk_size] for i in range(0, len(all_reads), chunk_size)]
    
    all_data = []
    
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_chunk = {
            executor.submit(
                process_read_chunk, 
                chunk, 
                primers_df, 
                window_size, 
                unaligned_only,
                max_distance,
                overlap_threshold=overlap_threshold  # Pass the parameter
            ): chunk for chunk in chunks
        }
        
        for future in tqdm(as_completed(future_to_chunk), total=len(chunks), desc="Processing chunks"):
            try:
                chunk_data = future.result()
                all_data.extend(chunk_data)
            except Exception as e:
                print(f"Error processing chunk: {e}")
    
    if not all_data:
        print("No data was processed successfully")
        return pd.DataFrame()
        
    return pd.DataFrame(all_data)

def process_read_chunk(chunk: List[Dict], primers_df: pd.DataFrame, 
                      window_size: int, unaligned_only: bool,
                      max_distance: int = 2, overlap_threshold: float = 0.8,
                      check_termini: bool = True, terminus_length: int = 10) -> List[Dict]:
    """
    Process a chunk of sequencing reads in parallel to identify primers.
    
    Args:
        ... (existing args)
        overlap_threshold: Minimum fraction of overlap required to consider primers as overlapping
    """
    chunk_data = []
    
    # Calculate max primer length once for efficiency
    max_primer_length = max(
        primers_df['Forward'].apply(len).max(),
        primers_df['Reverse'].apply(len).max()
    )
    
    # Adjusted window size including primer length
    effective_window = window_size + max_primer_length
    
    for read in chunk:
        # Handle different input types
        if isinstance(read, dict):  # Paired-end format
            read_sequence = read['sequence']
            read_name = read['name']
            is_unmapped = True  # Paired reads are already filtered
        else:  # Single-end format (pysam.AlignedSegment)
            if unaligned_only and not read.is_unmapped:
                continue
            if read.query_sequence is None:
                continue
            read_sequence = read.query_sequence
            read_name = read.query_name
            is_unmapped = read.is_unmapped

        if not read_sequence:
            continue
            
        read_length = len(read_sequence)
        
        # Get sequences from both ends of the read
        start_region = read_sequence[:min(effective_window, read_length)]
        end_start_pos = max(0, read_length - effective_window)
        end_region = read_sequence[end_start_pos:]
        
        # Process start primers
        start_results = find_primers_in_region(
            start_region, 
            primers_df, 
            window_size=window_size,
            max_distance=max_distance,
            check_termini=check_termini,
            terminus_length=terminus_length,
            overlap_threshold=overlap_threshold  # Pass the parameter
        )
        
        # Process end primers with reversed sequence
        end_results = find_primers_in_region(
            str(Seq(end_region).reverse_complement()),
            primers_df,
            window_size=window_size,
            max_distance=max_distance,
            check_termini=check_termini,
            terminus_length=terminus_length,
            overlap_threshold=overlap_threshold  # Pass the parameter
        )
        
        # Store results (rest of function remains the same)
        chunk_data.append({
            'Read_Name': read_name,
            'Start_Primers': ', '.join(start_results['full_matches']) if start_results['full_matches'] else 'None',
            'End_Primers': ', '.join(end_results['full_matches']) if end_results['full_matches'] else 'None',
            'Start_Terminal_Matches': ', '.join(start_results['terminal_matches']) if start_results['terminal_matches'] else 'None',
            'End_Terminal_Matches': ', '.join(end_results['terminal_matches']) if end_results['terminal_matches'] else 'None',
            'Read_Length': read_length,
            'Start_Region_Length': len(start_region),
            'End_Region_Length': len(end_region),
            'Is_Paired': isinstance(read, dict),
            'Is_Unmapped': is_unmapped
        })
    
    return chunk_data

def get_base_primer_name(primer_str):
    """
    Extract base primer name, handling both full matches and terminal matches.
    
    Args:
        primer_str (str): String containing primer match info
        
    Returns:
        str or None: Base primer name without orientation/terminal suffixes
    """
    if primer_str == 'None':
        return None
    primer_str = primer_str.split(',')[0].strip()
    # Handle terminal match format: Name_Orientation_Terminal_XXbp
    if '_Terminal_' in primer_str:
        return '_'.join(primer_str.split('_')[:-3])
    # Handle full match format: Name_Orientation
    return '_'.join(primer_str.split('_')[:-1])


def create_analysis_summary(result_df, primers_df, ignore_amplicon_size=False, debug=False, size_tolerance=0.10):
    """Create comprehensive summary of primer analysis results."""
    if result_df.empty:
        print("No reads to analyze in the results dataframe")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
    total_reads = len(result_df)
    
    # Existing helper functions remain the same
    def count_primers(primer_str):
        if primer_str == 'None':
            return 0
        return len([p for p in primer_str.split(',') if p.strip()])

    def get_terminal_length(match_str):
        if match_str == 'None':
            return 0
        try:
            # Extract length from format "PrimerName_Orientation_Terminal_XXbp"
            return int(match_str.split('_')[-1].replace('bp', ''))
        except:
            return 0
            
    # Add columns for analysis
    result_df['Start_Primer_Count'] = result_df['Start_Primers'].apply(count_primers)
    result_df['End_Primer_Count'] = result_df['End_Primers'].apply(count_primers)
    result_df['Start_Primer_Name'] = result_df['Start_Primers'].apply(get_base_primer_name)
    result_df['End_Primer_Name'] = result_df['End_Primers'].apply(get_base_primer_name)
    result_df['Start_Terminal_Count'] = result_df['Start_Terminal_Matches'].apply(count_primers)
    result_df['End_Terminal_Count'] = result_df['End_Terminal_Matches'].apply(count_primers)
    
    # Initialize columns
    result_df['Correct_Orientation'] = False
    result_df['Size_Compliant'] = False
    
    # Identify categories
    no_matches = result_df[
        (result_df['Start_Primers'] == 'None') & 
        (result_df['End_Primers'] == 'None') &
        (result_df['Start_Terminal_Matches'] == 'None') &
        (result_df['End_Terminal_Matches'] == 'None')
    ]
    
    single_terminal = result_df[
        (result_df['Start_Primers'] == 'None') & 
        (result_df['End_Primers'] == 'None') &
        ((result_df['Start_Terminal_Matches'] != 'None') & (result_df['End_Terminal_Matches'] == 'None') |
         (result_df['Start_Terminal_Matches'] == 'None') & (result_df['End_Terminal_Matches'] != 'None'))
    ]
    
    # Modified paired terminal analysis
    paired_terminal = result_df[
        (result_df['Start_Primers'] == 'None') & 
        (result_df['End_Primers'] == 'None') &
        (result_df['Start_Terminal_Matches'] != 'None') & 
        (result_df['End_Terminal_Matches'] != 'None')
    ]
    
    # Analyze length compliance for paired terminals
    paired_terminal_correct = pd.DataFrame()
    paired_terminal_incorrect = pd.DataFrame()
    
    if not paired_terminal.empty:
        for _, row in paired_terminal.iterrows():
            start_name = get_base_primer_name(row['Start_Terminal_Matches'])
            if start_name:
                expected_size = primers_df[primers_df['Name'] == start_name]['Size'].iloc[0]
                tolerance = expected_size * size_tolerance
                
                read_length = row['Read_Length']
                is_compliant = abs(read_length - expected_size) <= tolerance
                
                if is_compliant:
                    paired_terminal_correct = pd.concat([paired_terminal_correct, pd.DataFrame([row])])
                else:
                    paired_terminal_incorrect = pd.concat([paired_terminal_incorrect, pd.DataFrame([row])])
    
    hybrid_matches = result_df[
        ((result_df['Start_Primers'] != 'None') & (result_df['End_Primers'] == 'None') & (result_df['End_Terminal_Matches'] != 'None')) |
        ((result_df['Start_Primers'] == 'None') & (result_df['End_Primers'] != 'None') & (result_df['Start_Terminal_Matches'] != 'None'))
    ]
    
    single_end_only = result_df[
        ((result_df['Start_Primers'] != 'None') & (result_df['End_Primers'] == 'None') & (result_df['End_Terminal_Matches'] == 'None')) |
        ((result_df['Start_Primers'] == 'None') & (result_df['End_Primers'] != 'None') & (result_df['Start_Terminal_Matches'] == 'None'))
    ]
    
    matched_pairs = result_df[
        (result_df['Start_Primers'] != 'None') & 
        (result_df['End_Primers'] != 'None') &
        (result_df['Start_Primer_Count'] == 1) &
        (result_df['End_Primer_Count'] == 1) &
        (result_df['Start_Primer_Name'] == result_df['End_Primer_Name'])
    ].copy()
    
    multi_primer_pairs = result_df[
        (result_df['Start_Primers'] != 'None') & 
        (result_df['End_Primers'] != 'None') &
        ((result_df['Start_Primer_Count'] > 1) | (result_df['End_Primer_Count'] > 1))
    ]
    
    mismatched_pairs = result_df[
        (result_df['Start_Primers'] != 'None') & 
        (result_df['End_Primers'] != 'None') &
        (result_df['Start_Primer_Count'] == 1) &
        (result_df['End_Primer_Count'] == 1) &
        (result_df['Start_Primer_Name'] != result_df['End_Primer_Name'])
    ]
    
    if not matched_pairs.empty:
        matched_pairs['Correct_Orientation'] = matched_pairs.apply(
            lambda row: is_correct_orientation(row['Start_Primers'], row['End_Primers']), 
            axis=1
        )
        
        if not ignore_amplicon_size:
            matched_pairs['Size_Compliant'] = matched_pairs.apply(
                lambda row: is_size_compliant(row, primers_df, size_tolerance), 
                axis=1
            )
        else:
            matched_pairs['Size_Compliant'] = True

    # Create summary data
    summary_data = [
        {
            'Category': 'No primers or terminal matches detected',
            'Count': len(no_matches),
            'Percentage': (len(no_matches) / total_reads) * 100
        },
        {
            'Category': 'Single terminal match only',
            'Count': len(single_terminal),
            'Percentage': (len(single_terminal) / total_reads) * 100
        },
        {
            'Category': 'Paired terminal matches - correct length',
            'Count': len(paired_terminal_correct),
            'Percentage': (len(paired_terminal_correct) / total_reads) * 100
        },
        {
            'Category': 'Paired terminal matches - incorrect length',
            'Count': len(paired_terminal_incorrect),
            'Percentage': (len(paired_terminal_incorrect) / total_reads) * 100
        },
        {
            'Category': 'One full primer + one terminal match',
            'Count': len(hybrid_matches),
            'Percentage': (len(hybrid_matches) / total_reads) * 100
        },
        {
            'Category': 'Single-end primers only (no terminal match)',
            'Count': len(single_end_only),
            'Percentage': (len(single_end_only) / total_reads) * 100
        },
        {
            'Category': 'Multi-primer pairs (>1 primer at an end)',
            'Count': len(multi_primer_pairs),
            'Percentage': (len(multi_primer_pairs) / total_reads) * 100
        },
        {
            'Category': 'Mismatched primer pairs (different primers)',
            'Count': len(mismatched_pairs),
            'Percentage': (len(mismatched_pairs) / total_reads) * 100
        }
    ]
    
    if not matched_pairs.empty:
        correct_orient = matched_pairs[matched_pairs['Correct_Orientation']]
        if not ignore_amplicon_size:
            correct_orient_wrong_size = correct_orient[~correct_orient['Size_Compliant']]
            correct_orient_right_size = correct_orient[correct_orient['Size_Compliant']]
        else:
            correct_orient_wrong_size = pd.DataFrame()
            correct_orient_right_size = correct_orient

        summary_data.extend([
            {
                'Category': 'Matched pairs - correct orientation, wrong size',
                'Count': len(correct_orient_wrong_size),
                'Percentage': (len(correct_orient_wrong_size) / total_reads) * 100
            },
            {
                'Category': 'Matched pairs - correct orientation and size',
                'Count': len(correct_orient_right_size),
                'Percentage': (len(correct_orient_right_size) / total_reads) * 100
            }
        ])
    
    summary_df = pd.DataFrame(summary_data)
    summary_df['Percentage'] = summary_df['Percentage'].round(1)
    
    return summary_df, matched_pairs, mismatched_pairs

def is_correct_orientation(start_primers, end_primers):
    """Helper function to check if primers are in correct orientation"""
    if start_primers == 'None' or end_primers == 'None':
        return False
    
    start_orient = start_primers.split(',')[0].strip().split('_')[-1]
    end_orient = end_primers.split(',')[0].strip().split('_')[-1]
    
    return ((start_orient.startswith('Forward') and end_orient.startswith('Reverse')) or
            (start_orient.startswith('Reverse') and end_orient.startswith('Forward')))

def is_size_compliant(row, primers_df, size_tolerance=0.10):
    """
    Helper function to check if read length matches expected amplicon size
    """
    primer_name = get_base_primer_name(row['Start_Primers'])
    if primer_name is None:
        return False
        
    expected_size = primers_df[primers_df['Name'] == primer_name]['Size'].iloc[0]
    tolerance = expected_size * size_tolerance
    
    return abs(row['Read_Length'] - expected_size) <= tolerance

def is_illumina_data(bam_path: str) -> bool:
    """
    Check if the BAM file contains Illumina paired-end data by examining headers.
    
    Args:
        bam_path (str): Path to the BAM file
        
    Returns:
        bool: True if Illumina paired-end data is detected
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Check header for Illumina-specific tags
            header = bam.header
            
            # Look for Illumina platform in @PG or @RG tags
            pg_entries = header.get('PG', [])
            rg_entries = header.get('RG', [])
            
            for entry in pg_entries + rg_entries:
                platform = entry.get('PL', '').lower()
                if 'illumina' in platform:
                    return True
                    
            # Check first few reads for paired-end flags
            for i, read in enumerate(bam.fetch(until_eof=True)):
                if i >= 1000:  # Check first 1000 reads
                    break
                if read.is_paired:
                    return True
                    
            return False
            
    except Exception as e:
        print(f"Warning: Error checking for Illumina data: {e}")
        return False

def concatenate_paired_reads(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> Optional[str]:
    """
    Concatenate paired-end reads with 'N' separator.
    
    Args:
        read1: First read of the pair
        read2: Second read of the pair
        
    Returns:
        Optional[str]: Concatenated sequence or None if invalid
    """
    if not (read1 and read2 and read1.query_sequence and read2.query_sequence):
        return None
        
    # Add 'N' spacer between reads
    return read1.query_sequence + 'N' * 10 + read2.query_sequence

def process_paired_reads(bam_path: str, percentage: float, max_reads: int = 0) -> List[Dict]:
    """Process paired-end reads from BAM file with improved pair matching."""
    processed_reads = []
    pairs_dict = {}
    keep_probability = percentage / 100.0
    
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            reads_processed = 0
            pair_count = 0
            
            # First pass: collect read pairs
            print("Collecting read pairs...")
            for read in tqdm(bam.fetch(until_eof=True)):
                if not read.is_paired:
                    continue
                    
                reads_processed += 1
                
                # Apply downsampling at pair level
                if random.random() > keep_probability:
                    continue
                    
                if max_reads > 0 and pair_count >= max_reads:
                    break
                    
                qname = read.query_name
                
                if qname in pairs_dict:
                    pair = pairs_dict[qname]
                    # Make sure we have both reads and they're properly paired
                    if ((read.is_read1 and pair.is_read2) or 
                        (read.is_read2 and pair.is_read1)):
                        if read.is_read1:
                            read1, read2 = read, pair
                        else:
                            read1, read2 = pair, read
                            
                        if read1.query_sequence and read2.query_sequence:
                            concatenated_seq = concatenate_paired_reads(read1, read2)
                            if concatenated_seq:
                                processed_reads.append({
                                    'name': qname,
                                    'sequence': concatenated_seq,
                                    'is_paired': True
                                })
                                pair_count += 1
                    pairs_dict.pop(qname)  # Remove processed pair
                else:
                    pairs_dict[qname] = read
            
            print(f"Processed {pair_count} complete pairs out of {reads_processed} total reads")
            
    except Exception as e:
        print(f"Error processing paired reads: {e}")
        
    return processed_reads

def downsample_reads(bam_path: str, percentage: float, max_reads: int = 0) -> List[pysam.AlignedSegment]:
    """
    Downsample reads from a BAM file based on a percentage.
    
    Args:
        bam_path (str): Path to the BAM file
        percentage (float): Percentage of reads to keep (0.1-100.0)
        max_reads (int): Maximum number of reads to process (0 for all reads)
    
    Returns:
        List[pysam.AlignedSegment]: List of downsampled reads
    """

    if not (0.1 <= percentage <= 100.0):
        raise ValueError("Downsampling percentage must be between 0.1 and 100.0")
        
    # Check if data is Illumina paired-end
    is_paired = is_illumina_data(bam_path)
    
    if is_paired:
        print("Detected Illumina paired-end data - processing read pairs...")
        return process_paired_reads(bam_path, percentage, max_reads)
    else:
        print("Processing as single-end/long-read data...")
        try:
            bam_file = pysam.AlignmentFile(bam_path, "rb")
            
            # First pass to count total reads if needed
            total_reads = 0
            if max_reads == 0:
                print("Counting total reads...")
                for _ in tqdm(bam_file.fetch(until_eof=True)):
                    total_reads += 1
                bam_file.reset()
            else:
                total_reads = max_reads
                
            # Calculate number of reads to keep
            keep_probability = percentage / 100.0
            target_reads = int(total_reads * keep_probability)
            
            if max_reads > 0:
                target_reads = min(target_reads, max_reads)
                
            print(f"Targeting {target_reads} reads after {percentage}% downsampling")
            
            # Second pass to collect downsampled reads
            downsampled_reads = []
            reads_processed = 0
            
            for read in tqdm(bam_file.fetch(until_eof=True), total=total_reads):
                reads_processed += 1
                
                # Use reservoir sampling if we don't know total reads
                if max_reads == 0:
                    if len(downsampled_reads) < target_reads:
                        downsampled_reads.append(read)
                    else:
                        j = random.randint(0, reads_processed)
                        if j < target_reads:
                            downsampled_reads[j] = read
                # Otherwise use simple random sampling
                else:
                    if random.random() < keep_probability:
                        downsampled_reads.append(read)
                        if len(downsampled_reads) >= target_reads:
                            break
                
                if max_reads > 0 and reads_processed >= max_reads:
                    break
                    
            bam_file.close()
            
            print(f"Selected {len(downsampled_reads)} reads after downsampling")
            return downsampled_reads
            
        except Exception as e:
            print(f"Error during downsampling: {e}")
            return []

def parallel_analysis_pipeline(bam_path: str, primer_file: str, window_size: int = 20,
                             num_threads: int = 4, max_reads: int = 200, chunk_size: int = 50,
                             ignore_amplicon_size: bool = False,
                             max_distance: int = 2, 
                             downsample_percentage: float = 100.0,
                             unaligned_only: bool = False,
                             debug: bool = False,
                             size_tolerance: float = 0.10,
                             overlap_threshold: float = 0.8):
    """
    Complete analysis pipeline using parallel processing.
    
    Args:
        bam_path (str): Path to BAM file
        primer_file (str): Path to primer file
        window_size (int): Size of window to search for primers
        num_threads (int): Number of threads to use
        max_reads (int): Maximum number of reads to process
        chunk_size (int): Number of reads per chunk
        ignore_amplicon_size (bool): Whether to ignore amplicon size checks
        max_distance (int): Maximum Levenshtein distance for matching
        downsample_percentage (float): Percentage of reads to analyze
        unaligned_only (bool): Whether to process only unaligned reads
        debug (bool): Whether to print debug information
        size_tolerance (float): Size tolerance as fraction of expected size
        overlap_threshold (float): Minimum fraction of overlap for overlapping primers
    """
    print(f"Starting analysis with {num_threads} threads...")
    print(f"Using size tolerance of {size_tolerance:.1%}")
    
    try:
        # Get all results at once instead of processing in chunks
        result_df = bam_to_fasta_parallel(
            bam_path=bam_path,
            primer_file=primer_file,
            window_size=window_size,
            max_reads=max_reads,
            num_threads=num_threads,
            chunk_size=chunk_size,
            max_distance=max_distance,
            downsample_percentage=downsample_percentage,
            unaligned_only=unaligned_only,
            overlap_threshold=overlap_threshold
        )
        
        if result_df.empty:
            print("No results generated. Check input files and parameters.")
            return None
        
        print(f"\nProcessed {len(result_df)} reads successfully")
        
        primers_df, _ = load_primers(primer_file)
        
        # Process all data at once
        summary_df, matched_pairs, mismatched_pairs = create_analysis_summary(
            result_df.copy(),  # Use copy to prevent modifications
            primers_df,
            ignore_amplicon_size=ignore_amplicon_size,
            debug=debug,
            size_tolerance=size_tolerance
        )
        
        print("\nAnalysis Summary:")
        print(summary_df.to_string(index=False))
        
        return {
            'results': result_df.copy(),
            'summary': summary_df,
            'matched_pairs': matched_pairs.copy() if not matched_pairs.empty else pd.DataFrame(),
            'mismatched_pairs': mismatched_pairs.copy() if not mismatched_pairs.empty else pd.DataFrame()
        }
        
    except Exception as e:
        print(f"Error in analysis pipeline: {str(e)}")
        return None

def parse_arguments():
    """Parse command line arguments for URAdime."""
    parser = argparse.ArgumentParser(
        description="URAdime - Universal Read Analysis of DIMErs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input BAM file path"
    )
    
    parser.add_argument(
        "-p", "--primers",
        required=True,
        help="Tab-separated primer file containing columns: Name, Forward, Reverse, Size"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="uradime_results",
        help="Output prefix for result files"
    )
    
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads to use for parallel processing"
    )
    
    parser.add_argument(
        "-m", "--max-reads",
        type=int,
        default=0,
        help="Maximum number of reads to process (0 for all reads)"
    )
    
    parser.add_argument(
        "-c", "--chunk-size",
        type=int,
        default=50,
        help="Number of reads to process in each thread chunk"
    )
    
    parser.add_argument(
        "-u", "--unaligned-only",
        action="store_true",
        help="Process only unaligned reads"
    )
    
    parser.add_argument(
        "--max-distance",
        type=int,
        default=4,
        help="Maximum Levenshtein distance for primer matching"
    )
    
    parser.add_argument(
        "-w", "--window-size",
        type=int,
        default=30,
        help="Size of the window to search for primers at read ends"
    )

    parser.add_argument("--ignore-amplicon-size", 
        action="store_true", 
        help="Ignore amplicon size compliance checks"
        )

    parser.add_argument(
        "-d", "--downsample",
        type=float,
        default=100.0,
        help="Percentage of reads to randomly sample from the BAM file (0.1-100.0)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )

    parser.add_argument("--check-termini", 
        action="store_false", 
        help="Dont check for partial matches at read termini")
    
    parser.add_argument("--terminus-length", 
        type=int, 
        default=14, 
        help="Length of terminus to check for partial matches"
        )

    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=0.8,
        help="Minimum fraction of overlap required to consider primers as overlapping (0.0-1.0)"
    )

    parser.add_argument(
        "--size-tolerance",
        type=float,
        default=0.10,
        help="Size tolerance as fraction of expected amplicon size (e.g., 0.10 for 10%%)"
    )
        
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print detailed debug information for single-end primer reads"
    )
    
    return parser.parse_args()

def validate_inputs(args):
    """Validate input files and parameters."""
    if not os.path.exists(args.bam):
        raise FileNotFoundError(f"BAM file not found: {args.bam}")
    
    if not os.path.exists(args.primers):
        raise FileNotFoundError(f"Primer file not found: {args.primers}")

    if args.downsample <= 0 or args.downsample > 100:
        raise ValueError("Downsampling percentage must be between 0.1 and 100")
    
    try:
        primers_df = pd.read_csv(args.primers, sep="\t")
        required_columns = ['Name', 'Forward', 'Reverse', 'Size']
        missing_columns = [col for col in required_columns if col not in primers_df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns in primer file: {', '.join(missing_columns)}")
    except Exception as e:
        raise ValueError(f"Error reading primer file: {str(e)}")
    
    if args.max_distance < 0:
        raise ValueError("Maximum distance must be non-negative")
    
    if args.threads < 1:
        raise ValueError("Number of threads must be positive")
    
    if args.chunk_size < 1:
        raise ValueError("Chunk size must be positive")
        
    if args.window_size < 1:
        raise ValueError("Window size must be positive")
    
    if args.overlap_threshold < 0 or args.overlap_threshold > 1:
        raise ValueError("Overlap threshold must be between 0.0 and 1.0")

    if args.size_tolerance <= 0 or args.size_tolerance > 1:
        raise ValueError("Size tolerance must be between 0 and 1 (e.g., 0.10 for 10%)")

def create_primer_statistics(matched_pairs, primers_df, total_reads):
    """Create statistics for each primer pair."""
    if matched_pairs.empty:
        return pd.DataFrame()
        
    # Create a clean copy of the data
    matched_pairs = matched_pairs.copy()
    primers_df = primers_df.copy()
    
    # Ensure required columns exist
    if 'Start_Primer_Name' not in matched_pairs.columns:
        matched_pairs['Start_Primer_Name'] = matched_pairs['Start_Primers'].apply(get_base_primer_name)
        
    if 'Correct_Orientation' not in matched_pairs.columns:
        matched_pairs['Correct_Orientation'] = False
        
    if 'Size_Compliant' not in matched_pairs.columns:
        matched_pairs['Size_Compliant'] = False
    
    primer_stats = []
    
    for primer_name in primers_df['Name'].unique():
        # Get all reads where this primer appears
        primer_matches = matched_pairs[matched_pairs['Start_Primer_Name'] == primer_name]
        
        if len(primer_matches) == 0:
            continue
            
        total_appearances = len(primer_matches)
        stats = {
            'Primer_Name': primer_name,
            'Total_Appearances': total_appearances,
            'Percentage_of_Total_Reads': round((total_appearances / total_reads * 100), 2),
            'Correct_Orientation_Percentage': round(
                (primer_matches['Correct_Orientation'].astype(int).sum() / total_appearances * 100), 2
            ),
            'Size_Compliant_Percentage': round(
                (primer_matches['Size_Compliant'].astype(int).sum() / total_appearances * 100), 2
            ),
            'Correct_Orientation_and_Size_Percentage': round(
                (((primer_matches['Correct_Orientation'] == True) & 
                  (primer_matches['Size_Compliant'] == True)).sum() / total_appearances * 100), 2
            )
        }
        
        primer_stats.append(stats)
    
    return pd.DataFrame(primer_stats)
    
def save_results(results, output_prefix, primers_df):
    """Save analysis results to files with additional primer combination summaries."""
    os.makedirs(os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else '.', exist_ok=True)
    
    def create_primer_combination_summary(df, total_reads):
        """Helper function to create summary of primer combinations"""
        if df.empty:
            return pd.DataFrame()
            
        # Create a clean copy of the data
        df = df.copy()
        
        # Group by start and end primers to count occurrences
        summary = df.groupby(['Start_Primers', 'End_Primers'], observed=True).size().reset_index()
        summary.columns = ['Start_Primers', 'End_Primers', 'Occurrence_Count']
        
        # Calculate percentage of total reads
        summary['Percent_of_Total_Reads'] = (summary['Occurrence_Count'] / total_reads * 100).round(2)
        summary = summary.sort_values('Occurrence_Count', ascending=False)
        
        return summary
    
    # Get total reads analyzed
    total_reads = len(results['results'])
    
    # Save summary
    results['summary'].to_csv(f"{output_prefix}_summary.csv", index=False)
    
    # Save all results with primer details
    results['results'].to_csv(f"{output_prefix}_all_results.csv", index=False)
    
    # Save matched pairs and their summary
    if not results['matched_pairs'].empty:
        # Save full matched pairs data with all columns
        results['matched_pairs'].to_csv(f"{output_prefix}_matched_pairs.csv", index=False)
        
        # Create and save matched pairs summary
        matched_summary = create_primer_combination_summary(results['matched_pairs'], total_reads)
        if not matched_summary.empty:
            matched_summary.to_csv(f"{output_prefix}_matched_pairs_summary.csv", index=False)
        
        # Generate and save primer statistics
        primer_stats = create_primer_statistics(
            results['matched_pairs'],
            primers_df,
            total_reads
        )
        if not primer_stats.empty:
            primer_stats.to_csv(f"{output_prefix}_primer_statistics.csv", index=False)
    
    # Save mismatched pairs and their summary
    if not results['mismatched_pairs'].empty:
        # Save full mismatched pairs data
        results['mismatched_pairs'].to_csv(f"{output_prefix}_mismatched_pairs.csv", index=False)
        
        # Create and save mismatched pairs summary
        mismatched_summary = create_primer_combination_summary(results['mismatched_pairs'], total_reads)
        if not mismatched_summary.empty:
            mismatched_summary.to_csv(f"{output_prefix}_mismatched_pairs_summary.csv", index=False)
    
    # Save multi-primer pairs
    multi_primer_pairs = results['results'][
        (results['results']['Start_Primers'] != 'None') & 
        (results['results']['End_Primers'] != 'None') &
        ((results['results']['Start_Primers'].str.count(',') > 0) | 
         (results['results']['End_Primers'].str.count(',') > 0))
    ]
    
    if not multi_primer_pairs.empty:
        # Save full multi-primer pairs data
        multi_primer_pairs.to_csv(f"{output_prefix}_multi_primer_pairs.csv", index=False)
        
        # Create and save multi-primer pairs summary
        multi_primer_summary = create_primer_combination_summary(multi_primer_pairs, total_reads)
        if not multi_primer_summary.empty:
            multi_primer_summary.to_csv(f"{output_prefix}_multi_primer_pairs_summary.csv", index=False)
    
    # Save wrong size pairs and their summary
    wrong_size_pairs = results['results'][
        (results['results']['Start_Primers'] != 'None') & 
        (results['results']['End_Primers'] != 'None') &
        (results['results']['Start_Primer_Count'] == 1) &
        (results['results']['End_Primer_Count'] == 1) &
        (results['results']['Start_Primer_Name'] == results['results']['End_Primer_Name']) &
        ~results['results'].get('Size_Compliant', True)  # Handle cases where Size_Compliant might not exist
    ].copy()
    
    if not wrong_size_pairs.empty:
        # Save full wrong size pairs data
        wrong_size_pairs.to_csv(f"{output_prefix}_wrong_size_pairs.csv", index=False)
        
        # Create and save wrong size pairs summary
        wrong_size_summary = create_primer_combination_summary(wrong_size_pairs, total_reads)
        if not wrong_size_summary.empty:
            wrong_size_summary.to_csv(f"{output_prefix}_wrong_size_pairs_summary.csv", index=False)

def format_summary_table(df):
    """Format summary DataFrame as a pretty ASCII table."""
    # Get maximum lengths for each column
    cat_width = max(len(str(x)) for x in df['Category']) + 2
    count_width = max(len(str(x)) for x in df['Count']) + 2
    pct_width = max(len(f"{x:.1f}" if isinstance(x, float) else str(x)) for x in df['Percentage']) + 2

    # Create header
    header = (f"{'Category':<{cat_width}} {'Count':>{count_width}} {'Percentage':>{pct_width}}")
    separator = "=" * (cat_width + count_width + pct_width + 6)
    
    # Format each row
    rows = []
    for _, row in df.iterrows():
        rows.append(
            f"{str(row['Category']):<{cat_width}} "
            f"{str(row['Count']):>{count_width}} "
            f"{f'{float(row['Percentage']):.1f}':>{pct_width}}"
        )
    
    # Combine all parts
    table = f"\n{separator}\n{header}\n{separator}\n"
    table += "\n".join(rows)
    table += f"\n{separator}\n"
    
    return table

def main():
    """
    Main execution function for URAdime.
    
    This function coordinates the entire execution flow of the tool,
    from argument parsing to result output, handling errors and providing
    appropriate exit codes.
    
    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    print_banner()
    
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Validate inputs
        validate_inputs(args)
        
        if args.verbose:
            print("Starting URAdime analysis...")
            print(f"Input BAM: {args.bam}")
            print(f"Input primers: {args.primers}")
            print(f"Using {args.threads} threads")
            print(f"Window size: {args.window_size}")
            print(f"Max distance: {args.max_distance}")
            print(f"Size tolerance: {args.size_tolerance:.1%}")
            print(f"Overlap threshold: {args.overlap_threshold}")
            print(f"Downsampling to {args.downsample}% of reads")
        
        # Process BAM file and analyze
        results = parallel_analysis_pipeline(
            bam_path=args.bam,
            primer_file=args.primers,
            window_size=args.window_size,
            unaligned_only=args.unaligned_only,
            max_reads=args.max_reads,
            num_threads=args.threads,
            chunk_size=args.chunk_size,
            downsample_percentage=args.downsample,
            max_distance=args.max_distance,
            overlap_threshold=args.overlap_threshold,
            size_tolerance=args.size_tolerance,
            ignore_amplicon_size=args.ignore_amplicon_size,
            debug=args.debug
        )
        
        if results is None:
            print("No results generated. Check input files and parameters.")
            return 1
        
        if args.verbose:
            print(f"\nProcessed {len(results['results'])} reads successfully")
        
        # Load primers for saving results
        primers_df, _ = load_primers(args.primers)
        
        # Save results
        save_results(results, args.output, primers_df)
        
        # Print summary to console
        print("\nAnalysis Summary:")
        print(format_summary_table(results['summary']))
        
        if args.verbose:
            print(f"\nResults saved with prefix: {args.output}")
        
        return 0
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1