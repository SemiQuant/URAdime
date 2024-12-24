#!/usr/bin/env python3
import pysam
import random
from Bio.Seq import Seq
import pandas as pd
import subprocess
import os
from typing import List, Dict, Tuple
import click
from tabulate import tabulate
import numpy as np
from pathlib import Path

def create_random_primer(length: int = None, gc_content: float = None) -> str:
    """Create a random primer with specified length and GC content"""
    if length is None:
        length = random.randint(18, 25)  # Common primer lengths
    if gc_content is None:
        gc_content = random.uniform(0.4, 0.6)  # Common GC content range
    
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(['G', 'C']))
        else:
            bases.append(random.choice(['A', 'T']))
    return ''.join(bases)

def create_test_primers(num_pairs: int = None, 
                       min_size: int = None, 
                       max_size: int = None) -> pd.DataFrame:
    """Create test primer pairs with varying parameters"""
    if num_pairs is None:
        num_pairs = random.randint(3, 10)
    if min_size is None:
        min_size = 100
    if max_size is None:
        max_size = 500
    
    primers = {
        'Name': [f'Pair{i+1}' for i in range(num_pairs)],
        'Forward': [create_random_primer() for _ in range(num_pairs)],
        'Reverse': [create_random_primer() for _ in range(num_pairs)],
        'Size': [random.randint(min_size, max_size) for _ in range(num_pairs)]
    }
    return pd.DataFrame(primers)

def generate_read(length: int = 100, gc_content: float = 0.5) -> str:
    """Generate random DNA sequence with specified GC content"""
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(['G', 'C']))
        else:
            bases.append(random.choice(['A', 'T']))
    return ''.join(bases)

def get_valid_size(size: int, primers_df: pd.DataFrame, row_idx: int) -> int:
    """Get size within tolerance (±10%)"""
    # Get primer lengths for this pair
    row = primers_df.iloc[row_idx]
    fwd_len = len(row['Forward'])
    rev_len = len(row['Reverse'])
    
    # The target size in primers.tsv already includes primer lengths
    # So we need to subtract primer lengths to get the insert size
    insert_size = size - (fwd_len + rev_len)
    
    # Generate a size within ±9.5% to be safe within URAdime's 10% tolerance
    min_size = int(insert_size * 0.905)
    max_size = int(insert_size * 1.095)
    
    # Generate a size that will result in the total read length being within tolerance
    valid_size = random.randint(min_size, max_size)
    total_size = valid_size + fwd_len + rev_len
    
    # Verify that the total size is within tolerance of the target size
    while abs(total_size - size) / size > 0.095:
        valid_size = random.randint(min_size, max_size)
        total_size = valid_size + fwd_len + rev_len
    
    return valid_size

def get_invalid_size(size: int, primers_df: pd.DataFrame, row_idx: int) -> int:
    """Get size outside tolerance (beyond ±10%)"""
    # Get primer lengths for this pair
    row = primers_df.iloc[row_idx]
    fwd_len = len(row['Forward'])
    rev_len = len(row['Reverse'])
    
    # Calculate insert size without primers
    insert_size = size - (fwd_len + rev_len)
    
    if random.random() < 0.5:
        # Shorter than expected: less than 90% of target
        invalid_size = max(20, int(insert_size * random.uniform(0.5, 0.89)))
        total_size = invalid_size + fwd_len + rev_len
        
        # Verify that the total size is outside tolerance
        while abs(total_size - size) / size <= 0.105:
            invalid_size = max(20, int(insert_size * random.uniform(0.5, 0.89)))
            total_size = invalid_size + fwd_len + rev_len
    else:
        # Longer than expected: more than 110% of target
        invalid_size = int(insert_size * random.uniform(1.11, 1.5))
        total_size = invalid_size + fwd_len + rev_len
        
        # Verify that the total size is outside tolerance
        while abs(total_size - size) / size <= 0.105:
            invalid_size = int(insert_size * random.uniform(1.11, 1.5))
            total_size = invalid_size + fwd_len + rev_len
    
    return invalid_size

def create_dataset_params() -> Dict:
    """Create random parameters for a dataset"""
    return {
        'num_primers': random.randint(5, 10),
        'min_size': random.randint(100, 200),
        'max_size': random.randint(300, 1000),
        'total_reads': random.randint(800, 1200),
        'gc_content': random.uniform(0.4, 0.6),
        'category_ratios': {
            'no_primers': random.uniform(0.06, 0.07),
            'correct_pairs': random.uniform(0.29, 0.31),  # ~30% for matched pairs
            'wrong_size': random.uniform(0.06, 0.07),
            'single_end': random.uniform(0.06, 0.07),
            'mismatched_pairs': random.uniform(0.06, 0.07),
            'full_terminal_match_correct': random.uniform(0.06, 0.07),
            'full_terminal_match_wrong': random.uniform(0.06, 0.07),
            'single_terminal_correct': random.uniform(0.06, 0.07),
            'single_terminal_wrong': random.uniform(0.06, 0.07),
            'paired_terminal_correct': random.uniform(0.06, 0.07),
            'paired_terminal_wrong': random.uniform(0.06, 0.07),
            'multi_primer': random.uniform(0.06, 0.07)
        }
    }

def generate_non_matching_sequence(length: int, avoid_seq: str, gc_content: float = 0.5) -> str:
    """Generate a sequence that's guaranteed to be different from any primer sequence"""
    # Create a sequence that's completely different from the primer
    # by using the complement of each base
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement_map[base] for base in avoid_seq[:length])

def get_valid_size_with_terminal(size: int, primers_df: pd.DataFrame, index: int) -> int:
    """Get a valid size that accounts for terminal matches"""
    row = primers_df.iloc[index]
    # For terminal matches with one full primer, we need to account for:
    # - One full primer length
    # - One 15bp terminal match
    full_primer_len = len(row['Forward'])  # Using forward as the full primer
    terminal_len = 15  # Standard terminal match length
    
    # Calculate the target insert size (excluding primers)
    target_insert = size - (full_primer_len + terminal_len)
    
    # Calculate tolerance range (9.5% to be safe within URAdime's 10% tolerance)
    min_size = int(target_insert * 0.905)
    max_size = int(target_insert * 1.095)
    
    # Generate a size that will result in the total read length being within tolerance
    valid_size = random.randint(min_size, max_size)
    total_size = valid_size + full_primer_len + terminal_len
    
    # Verify that the total size is within tolerance of the target size
    while abs(total_size - size) / size > 0.095:
        valid_size = random.randint(min_size, max_size)
        total_size = valid_size + full_primer_len + terminal_len
    
    return valid_size

def get_invalid_size_with_terminal(size: int, primers_df: pd.DataFrame, index: int) -> int:
    """Get an invalid size that accounts for terminal matches"""
    row = primers_df.iloc[index]
    # For terminal matches with one full primer, we need to account for:
    # - One full primer length
    # - One 15bp terminal match
    full_primer_len = len(row['Forward'])  # Using forward as the full primer
    terminal_len = 15  # Standard terminal match length
    
    # Calculate the target insert size (excluding primers)
    target_insert = size - (full_primer_len + terminal_len)
    
    # Generate a size that's either too small or too large
    if random.random() < 0.5:
        # Too small: less than 30% of target
        min_size = max(30, int(target_insert * 0.2))  # At least 30 bases or 20% of target size
        max_size = int(target_insert * 0.29)  # Just under 30%
    else:
        # Too large: more than 110% of target
        min_size = int(target_insert * 1.2)
        max_size = int(target_insert * 3.0)
    
    # Generate an invalid size
    invalid_size = random.randint(min_size, max_size)
    total_size = invalid_size + full_primer_len + terminal_len
    
    # Verify that the total size is outside tolerance
    while abs(total_size - size) / size <= 0.105:
        if random.random() < 0.5:
            invalid_size = max(30, int(target_insert * random.uniform(0.2, 0.29)))
        else:
            invalid_size = int(target_insert * random.uniform(1.2, 3.0))
        total_size = invalid_size + full_primer_len + terminal_len
    
    return invalid_size

def create_test_reads(primers_df: pd.DataFrame, params: Dict) -> List[Dict]:
    """Create test reads for each category with varying parameters"""
    reads = []
    total_reads = params['total_reads']
    gc_content = params['gc_content']
    
    # Calculate counts based on ratios
    category_counts = {
        cat: max(1, int(ratio * total_reads))
        for cat, ratio in params['category_ratios'].items()
    }
    
    # Adjust to ensure total matches
    total_count = sum(category_counts.values())
    if total_count != total_reads:
        # Adjust the largest category to make up the difference
        largest_cat = max(category_counts.items(), key=lambda x: x[1])[0]
        category_counts[largest_cat] += (total_reads - total_count)
    
    categories = {
        'no_primers': {
            'count': category_counts['no_primers'],
            'category': '🟥 No primers or terminal matches detected',
            'generator': lambda i, row: generate_read(
                random.randint(params['min_size'], params['max_size']),
                gc_content=gc_content
            )
        },
        'correct_pairs': {
            'count': category_counts['correct_pairs'],
            'category': '🟩 Matched pairs - correct orientation and size',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_valid_size(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse']).reverse_complement())
            )
        },
        'wrong_size': {
            'count': category_counts['wrong_size'],
            'category': '🟧 Matched pairs - correct orientation, wrong size',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_invalid_size(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse']).reverse_complement())
            )
        },
        'single_end': {
            'count': category_counts['single_end'],
            'category': '🟨 Single-end primers only (no terminal match)',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_valid_size(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                )
            )
        },
        'mismatched_pairs': {
            'count': category_counts['mismatched_pairs'],
            'category': '🟥 Mismatched primer pairs (different primers)',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_valid_size(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(primers_df.iloc[(i + 1) % len(primers_df)]['Reverse']).reverse_complement())
            )
        },
        'full_terminal_match_correct': {
            'count': category_counts['full_terminal_match_correct'],
            'category': '🟨 One full primer + one terminal match - correct size',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_valid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'full_terminal_match_wrong': {
            'count': category_counts['full_terminal_match_wrong'],
            'category': '🟥 One full primer + one terminal match - wrong size',
            'generator': lambda i, row: (
                row['Forward'] +
                generate_read(
                    get_invalid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'single_terminal_correct': {
            'count': category_counts['single_terminal_correct'],
            'category': '🟨 Single terminal match only - correct size',
            'generator': lambda i, row: (
                generate_read(
                    get_valid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'single_terminal_wrong': {
            'count': category_counts['single_terminal_wrong'],
            'category': '🟥 Single terminal match only - wrong size',
            'generator': lambda i, row: (
                generate_read(
                    get_invalid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'paired_terminal_correct': {
            'count': category_counts['paired_terminal_correct'],
            'category': '🟨 Paired terminal matches - correct size',
            'generator': lambda i, row: (
                row['Forward'][:15] +  # First 15 bases of forward primer
                generate_read(
                    get_valid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'paired_terminal_wrong': {
            'count': category_counts['paired_terminal_wrong'],
            'category': '🟥 Paired terminal matches - wrong size',
            'generator': lambda i, row: (
                row['Forward'][:15] +  # First 15 bases of forward primer
                generate_read(
                    get_invalid_size_with_terminal(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse'][:15]).reverse_complement())  # First 15 bases from 5' end
            )
        },
        'multi_primer': {
            'count': category_counts['multi_primer'],
            'category': '🟥 Multi-primer pairs (>1 primer at an end)',
            'generator': lambda i, row: (
                row['Forward'] +
                primers_df.iloc[(i + 1) % len(primers_df)]['Forward'] +
                generate_read(
                    get_valid_size(row['Size'], primers_df, i % len(primers_df)),
                    gc_content=gc_content
                ) +
                str(Seq(row['Reverse']).reverse_complement())
            )
        }
    }
    
    # Generate reads for each category
    for cat_name, cat_info in categories.items():
        for i in range(cat_info['count']):
            row = primers_df.iloc[i % len(primers_df)]
            try:
                sequence = cat_info['generator'](i, row)
                reads.append({
                    'name': f'{cat_name}_{row["Name"]}_{i}',
                    'sequence': sequence,
                    'category': cat_info['category']
                })
            except Exception as e:
                print(f"Error generating read for {cat_name}: {e}")
    
    # Shuffle reads to avoid any potential bias
    random.shuffle(reads)
    return reads[:total_reads]  # Ensure we return exactly the requested number of reads

def create_bam_file(reads: List[Dict], output_bam: str):
    """Create BAM file from generated reads"""
    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1000, 'SN': 'test_ref'}]}
    
    with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
        for read in reads:
            a = pysam.AlignedSegment()
            a.query_name = read['name']
            a.query_sequence = read['sequence']
            a.flag = 4  # unmapped
            a.reference_id = -1
            a.reference_start = 0
            a.mapping_quality = 0
            a.query_qualities = pysam.qualitystring_to_array("I" * len(read['sequence']))
            outf.write(a)
    
    # Index the BAM file
    pysam.index(output_bam)

def format_comparison_table(comparison_df: pd.DataFrame) -> str:
    """Format comparison results as a pretty table"""
    # Sort by category type (emoji) and then by name
    comparison_df['sort_key'] = comparison_df['Category'].str[0]
    comparison_df = comparison_df.sort_values('sort_key')
    
    # Add agreement column
    comparison_df['Agreement'] = (comparison_df['Count_expected'] == comparison_df['Count_actual'])
    comparison_df['Agreement'] = comparison_df['Agreement'].map({True: '✅', False: '❌'})
    
    # Calculate percentage difference
    comparison_df['Diff'] = comparison_df['Count_actual'] - comparison_df['Count_expected']
    comparison_df['Diff'] = comparison_df['Diff'].astype(int)
    
    # Format table with clear sections
    table = "\n=== Test Results ===\n\n"
    
    # Group by emoji and track section stats
    section_stats = {}
    for emoji in ['🟩', '🟧', '🟨', '🟥']:
        emoji_group = comparison_df[comparison_df['Category'].str.startswith(emoji)]
        if not emoji_group.empty:
            table_data = []
            section_matches = 0
            for _, row in emoji_group.iterrows():
                table_data.append([
                    row['Category'],
                    int(row['Count_expected']),
                    int(row['Count_actual']),
                    f"{int(row['Diff']):+d}",
                    f"{row['Percentage_expected']:.1f}%",
                    f"{row['Percentage_actual']:.1f}%",
                    row['Agreement']
                ])
                if row['Agreement'] == '✅':
                    section_matches += 1
            
            # Calculate section agreement rate
            section_total = len(emoji_group)
            section_agreement = (section_matches / section_total * 100) if section_total > 0 else 100
            section_stats[emoji] = {
                'matches': section_matches,
                'total': section_total,
                'rate': section_agreement
            }
            
            # Add section header with agreement rate
            section_name = {
                '🟩': 'Green',
                '🟧': 'Orange',
                '🟨': 'Yellow',
                '🟥': 'Red'
            }[emoji]
            
            table += f"=== {section_name} Section: {section_agreement:.1f}% Agreement ({section_matches}/{section_total}) ===\n"
            
            # Add section to table
            section = tabulate(
                table_data,
                headers=['Category', 'Exp', 'Act', 'Diff', 'Exp %', 'Act %', 'Match'],
                tablefmt='simple',
                showindex=False
            )
            table += f"{section}\n\n"
    
    # Calculate and add overall summary
    total_matches = sum(stats['matches'] for stats in section_stats.values())
    total_categories = sum(stats['total'] for stats in section_stats.values())
    agreement_rate = (total_matches / total_categories) * 100
    
    table += "=== Overall Summary ===\n"
    table += f"Total Agreement Rate: {agreement_rate:.1f}% ({total_matches}/{total_categories} categories)\n"
    
    # Add section-specific summary
    table += "\nSection-specific Agreement Rates:\n"
    for emoji, name in [('🟩', 'Green'), ('🟧', 'Orange'), ('🟨', 'Yellow'), ('🟥', 'Red')]:
        stats = section_stats[emoji]
        table += f"{name}: {stats['rate']:.1f}% ({stats['matches']}/{stats['total']} categories)\n"
    
    return table

def save_dataset_params(params: Dict, output_dir: str, dataset_num: int):
    """Save dataset parameters to a JSON file"""
    import json
    params_file = os.path.join(output_dir, f'dataset_{dataset_num}_params.json')
    with open(params_file, 'w') as f:
        json.dump(params, f, indent=2)

def create_final_tally(output_dir: str) -> Tuple[str, pd.DataFrame]:
    """Create a final tally table summarizing all datasets"""
    all_datasets = []
    
    # Collect all dataset results
    for dataset_dir in Path(output_dir).glob('dataset_*'):
        comparison_file = dataset_dir / 'comparison.csv'
        if comparison_file.exists():
            df = pd.read_csv(comparison_file)
            df['Dataset'] = dataset_dir.name
            all_datasets.append(df)
    
    if not all_datasets:
        return "No datasets found!", pd.DataFrame()
    
    # Combine all results
    combined_df = pd.concat(all_datasets)
    
    # Calculate statistics per category
    stats = []
    for category in combined_df['Category'].unique():
        category_data = combined_df[combined_df['Category'] == category]
        total_expected = category_data['Count_expected'].sum()
        total_actual = category_data['Count_actual'].sum()
        num_datasets = len(category_data)
        perfect_matches = sum(category_data['Count_expected'] == category_data['Count_actual'])
        
        stats.append({
            'Category': category,
            'Total Expected': total_expected,
            'Total Actual': total_actual,
            'Difference': total_actual - total_expected,
            'Perfect Matches': f"{perfect_matches}/{num_datasets}",
            'Match Rate': f"{(total_actual/total_expected*100):.1f}%" if total_expected > 0 else "0.0%"
        })
    
    # Create summary table
    stats_df = pd.DataFrame(stats)
    stats_df['sort_key'] = stats_df['Category'].str[0]
    stats_df = stats_df.sort_values('sort_key').drop('sort_key', axis=1)
    
    table = "\n=== Final Tally Across All Datasets ===\n\n"
    
    # Add section summaries
    section_totals = {
        '🟩': {'expected': 0, 'actual': 0, 'perfect': 0, 'total': 0},
        '🟧': {'expected': 0, 'actual': 0, 'perfect': 0, 'total': 0},
        '🟨': {'expected': 0, 'actual': 0, 'perfect': 0, 'total': 0},
        '🟥': {'expected': 0, 'actual': 0, 'perfect': 0, 'total': 0}
    }
    
    section_stats = []
    
    for emoji in ['🟩', '🟧', '🟨', '🟥']:
        section_data = stats_df[stats_df['Category'].str.startswith(emoji)]
        if not section_data.empty:
            section_name = {
                '🟩': 'Green',
                '🟧': 'Orange',
                '🟨': 'Yellow',
                '🟥': 'Red'
            }[emoji]
            
            # Update section totals
            section_totals[emoji]['expected'] = section_data['Total Expected'].sum()
            section_totals[emoji]['actual'] = section_data['Total Actual'].sum()
            section_totals[emoji]['perfect'] = sum([int(x.split('/')[0]) for x in section_data['Perfect Matches']])
            section_totals[emoji]['total'] = sum([int(x.split('/')[1]) for x in section_data['Perfect Matches']])
            
            table += f"=== {section_name} Section ===\n"
            section_table = tabulate(
                section_data,
                headers='keys',
                tablefmt='simple',
                showindex=False
            )
            table += f"{section_table}\n\n"
            
            # Add section data to DataFrame
            for _, row in section_data.iterrows():
                section_stats.append(row.to_dict())
    
    # Add overall statistics in the same format
    table += "=== Overall Statistics ===\n"
    overall_stats = []
    
    # Add section totals
    for emoji, name in [('🟩', 'Green Total'), ('🟧', 'Orange Total'), ('🟨', 'Yellow Total'), ('🟥', 'Red Total')]:
        totals = section_totals[emoji]
        if totals['total'] > 0:  # Only add if section has data
            stats = {
                'Category': name,
                'Total Expected': totals['expected'],
                'Total Actual': totals['actual'],
                'Difference': totals['actual'] - totals['expected'],
                'Perfect Matches': f"{totals['perfect']}/{totals['total']}",
                'Match Rate': f"{(totals['actual']/totals['expected']*100):.1f}%" if totals['expected'] > 0 else "0.0%"
            }
            overall_stats.append(stats)
            section_stats.append(stats)
    
    # Add grand total with updated match rate calculation
    grand_total = {
        'Category': 'Grand Total',
        'Total Expected': sum(x['expected'] for x in section_totals.values()),
        'Total Actual': sum(x['actual'] for x in section_totals.values()),
        'Difference': sum(x['actual'] - x['expected'] for x in section_totals.values()),
        'Perfect Matches': f"{sum(x['perfect'] for x in section_totals.values())}/{sum(x['total'] for x in section_totals.values())}",
        'Match Rate': f"{(sum(x['actual'] for x in section_totals.values())/sum(x['expected'] for x in section_totals.values())*100):.1f}%"
    }
    overall_stats.append(grand_total)
    section_stats.append(grand_total)
    
    overall_table = tabulate(
        overall_stats,
        headers='keys',
        tablefmt='simple',
        showindex=False
    )
    table += f"{overall_table}\n"
    
    return table, pd.DataFrame(section_stats)

def create_dataset_summary(output_dir: str) -> str:
    """Create a summary table of dataset parameters"""
    import json
    
    summaries = []
    for dataset_dir in Path(output_dir).glob('dataset_*'):
        params_file = dataset_dir / f'{dataset_dir.name}_params.json'
        primers_file = dataset_dir / 'primers.tsv'
        
        if params_file.exists() and primers_file.exists():
            # Load parameters
            with open(params_file) as f:
                params = json.load(f)
            
            # Load primers
            primers_df = pd.read_csv(primers_file, sep='\t')
            
            # Calculate average GC content of primers
            gc_content = []
            for _, row in primers_df.iterrows():
                fwd_gc = (row['Forward'].count('G') + row['Forward'].count('C')) / len(row['Forward'])
                rev_gc = (row['Reverse'].count('G') + row['Reverse'].count('C')) / len(row['Reverse'])
                gc_content.extend([fwd_gc, rev_gc])
            avg_primer_gc = np.mean(gc_content)
            
            # Calculate average primer length
            avg_primer_len = np.mean([len(row['Forward']) + len(row['Reverse']) 
                                    for _, row in primers_df.iterrows()])
            
            # Get the category ratio (use the first category in the ratios dict)
            category = next(iter(params['category_ratios']))
            ratio = params['category_ratios'][category]
            
            summaries.append({
                'Dataset': dataset_dir.name,
                'Primer Sets': len(primers_df),
                'Min Size': params['min_size'],
                'Max Size': params['max_size'],
                'Total Reads': params['total_reads'],
                'Read GC Content': f"{params['gc_content']:.2f}",
                'Avg Primer GC': f"{avg_primer_gc:.2f}",
                'Avg Primer Length': f"{avg_primer_len:.1f}",
                'Category': category,
                'Category Ratio %': f"{ratio*100:.1f}%"
            })
    
    if not summaries:
        return "No datasets found!"
    
    # Create summary table
    summary_df = pd.DataFrame(summaries)
    table = "\n=== Dataset Parameters Summary ===\n\n"
    table += tabulate(
        summary_df,
        headers='keys',
        tablefmt='simple',
        showindex=False
    )
    
    return table, summary_df

def save_to_excel(output_dir: str, tally_df: pd.DataFrame, summary_df: pd.DataFrame):
    """Save tally and summary tables to an Excel file"""
    excel_path = os.path.join(output_dir, 'test_datasets_report.xlsx')
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        # Write the tally table
        tally_df.to_excel(writer, sheet_name='Category Tally', index=False)
        
        # Write the summary table
        summary_df.to_excel(writer, sheet_name='Dataset Summary', index=False)
        
        # Auto-adjust column widths
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            for column in worksheet.columns:
                max_length = 0
                column = [cell for cell in column]
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = (max_length + 2)
                worksheet.column_dimensions[column[0].column_letter].width = adjusted_width

def run_uradime_analysis(dataset_dir):
    """Run URAdime analysis on the generated dataset."""
    bam_file = os.path.join(dataset_dir, "reads.bam")
    primers_file = os.path.join(dataset_dir, "primers.csv")
    
    cmd = [
        "uradime",
        "--input", bam_file,
        "--primers", primers_file,
        "--max-distance", "0",
        "--threads", "8"
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running URAdime: {e}")
        return None
        
    results_file = "uradime_results_summary.csv"
    if not os.path.exists(results_file):
        print(f"Results file not found: {results_file}")
        return None
        
    return pd.read_csv(results_file)

@click.command()
@click.option('--output-dir', default='test_data', help='Output directory for test files')
@click.option('--num-datasets', default=10, help='Number of datasets to generate')
@click.option('--verbose', is_flag=True, help='Show detailed output')
def main(output_dir: str, num_datasets: int, verbose: bool):
    """Generate multiple test datasets with varying parameters"""
    # Clear the output directory if it exists
    if os.path.exists(output_dir):
        import shutil
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    for dataset_num in range(1, num_datasets + 1):
        dataset_dir = os.path.join(output_dir, f'dataset_{dataset_num}')
        os.makedirs(dataset_dir, exist_ok=True)
        
        click.echo(f"\nGenerating dataset {dataset_num}/{num_datasets}")
        with click.progressbar(length=4, label='Generating test data') as bar:
            # Generate random parameters for this dataset
            params = create_dataset_params()
            save_dataset_params(params, dataset_dir, dataset_num)
            
            # Create primers with varying parameters
            primers_df = create_test_primers(
                num_pairs=params['num_primers'],
                min_size=params['min_size'],
                max_size=params['max_size']
            )
            primers_df.to_csv(os.path.join(dataset_dir, 'primers.tsv'), sep='\t', index=False)
            bar.update(1)
            
            # Generate reads
            reads = create_test_reads(primers_df, params)
            
            # Debug output for problematic categories
            if verbose:
                click.echo("\nSample reads from problematic categories:")
                problem_cats = ['full_terminal_match_correct', 'full_terminal_match_wrong', 
                              'single_terminal_correct', 'single_terminal_wrong',
                              'paired_terminal_correct', 'paired_terminal_wrong',
                              'multi_primer']
                for cat in problem_cats:
                    cat_reads = [r for r in reads if r['name'].startswith(cat)]
                    if cat_reads:
                        read = cat_reads[0]
                        click.echo(f"\n{read['category']}:")
                        click.echo(f"Name: {read['name']}")
                        click.echo(f"Sequence length: {len(read['sequence'])}")
                        click.echo(f"First 20 bases: {read['sequence'][:20]}...")
                        click.echo(f"Last 20 bases: ...{read['sequence'][-20:]}")
            
            bam_path = os.path.join(dataset_dir, 'reads.bam')
            create_bam_file(reads, bam_path)
            bar.update(1)
            
            # Create expected results
            categories = set(r['category'] for r in reads)
            expected_results = []
            total_reads = len(reads)
            
            for cat in categories:
                count = len([r for r in reads if r['category'] == cat])
                expected_results.append({
                    'Category': cat,
                    'Count': count,
                    'Percentage': (count / total_reads * 100)
                })
            
            expected_df = pd.DataFrame(expected_results)
            expected_df.to_csv(os.path.join(dataset_dir, 'expected_results.csv'), index=False)
            bar.update(1)
            
            # Run URAdime
            click.echo("Running URAdime analysis...")
            cmd = [
                'uradime',
                '-b', bam_path,
                '-p', os.path.join(dataset_dir, 'primers.tsv'),
                '-o', os.path.join(dataset_dir, 'uradime_results'),
                '--max-distance', '0'
            ]
            
            subprocess.run(cmd, capture_output=not verbose)
            bar.update(1)
            
            # Load actual results and compare
            actual_df = pd.read_csv(os.path.join(dataset_dir, 'uradime_results_summary.csv'))
            
            # Merge and display comparison
            comparison = pd.merge(
                expected_df, 
                actual_df, 
                on='Category', 
                suffixes=('_expected', '_actual'),
                how='outer'
            ).fillna(0)
            
            # Save comparison
            comparison.to_csv(os.path.join(dataset_dir, 'comparison.csv'), index=False)
            
            if verbose:
                click.echo(f"\nDataset {dataset_num} Results:")
                click.echo(format_comparison_table(comparison))
    
    # Create and display final tally
    click.echo("\nFinal Tally:")
    tally_table, tally_df = create_final_tally(output_dir)
    click.echo(tally_table)
    
    # Create and display dataset summary
    click.echo("\nDataset Summary:")
    summary_table, summary_df = create_dataset_summary(output_dir)
    # click.echo(summary_table)
    
    # Save both tables to Excel
    save_to_excel(output_dir, tally_df, summary_df)
    click.echo(f"\nTables have been saved to {os.path.join(output_dir, 'test_datasets_report.xlsx')}")

if __name__ == '__main__':
    main() 