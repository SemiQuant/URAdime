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
        'num_primers': random.randint(3, 10),
        'min_size': random.randint(100, 200),
        'max_size': random.randint(300, 1000),
        'total_reads': random.randint(800, 1200),
        'gc_content': random.uniform(0.4, 0.6),
        'category_ratios': {
            'no_primers': random.uniform(0.3, 0.4),
            'correct_pairs': random.uniform(0.01, 0.05),
            'wrong_size': random.uniform(0.4, 0.5),
            'single_end': random.uniform(0.1, 0.2),
            'mismatched_pairs': random.uniform(0.001, 0.01)
        }
    }

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
                random.randint(40, 60),
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
                    get_invalid_size(row['Size'], primers_df, i % len(primers_df)),
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
    return reads

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

def create_final_tally(output_dir: str) -> str:
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
        return "No datasets found!"
    
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
            'Match Rate': f"{perfect_matches/num_datasets*100:.1f}%"
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
    
    # Add overall statistics in the same format
    table += "=== Overall Statistics ===\n"
    overall_stats = []
    
    # Add section totals
    for emoji, name in [('🟩', 'Green Total'), ('🟧', 'Orange Total'), ('🟨', 'Yellow Total'), ('🟥', 'Red Total')]:
        totals = section_totals[emoji]
        if totals['total'] > 0:  # Only add if section has data
            overall_stats.append({
                'Category': name,
                'Total Expected': totals['expected'],
                'Total Actual': totals['actual'],
                'Difference': totals['actual'] - totals['expected'],
                'Perfect Matches': f"{totals['perfect']}/{totals['total']}",
                'Match Rate': f"{(totals['perfect']/totals['total']*100):.1f}%" if totals['total'] > 0 else "0.0%"
            })
    
    # Add grand total
    grand_total = {
        'Category': 'Grand Total',
        'Total Expected': sum(x['expected'] for x in section_totals.values()),
        'Total Actual': sum(x['actual'] for x in section_totals.values()),
        'Difference': sum(x['actual'] - x['expected'] for x in section_totals.values()),
        'Perfect Matches': f"{sum(x['perfect'] for x in section_totals.values())}/{sum(x['total'] for x in section_totals.values())}",
        'Match Rate': f"{(sum(x['perfect'] for x in section_totals.values())/sum(x['total'] for x in section_totals.values())*100):.1f}%"
    }
    overall_stats.append(grand_total)
    
    overall_table = tabulate(
        overall_stats,
        headers='keys',
        tablefmt='simple',
        showindex=False
    )
    table += f"{overall_table}\n"
    
    return table

@click.command()
@click.option('--output-dir', default='test_data', help='Output directory for test files')
@click.option('--num-datasets', default=10, help='Number of datasets to generate')
@click.option('--verbose', is_flag=True, help='Show detailed output')
def main(output_dir: str, num_datasets: int, verbose: bool):
    """Generate multiple test datasets with varying parameters"""
    os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Display final tally
    click.echo("\nFinal Tally:")
    click.echo(create_final_tally(output_dir))

if __name__ == '__main__':
    main() 