#!/usr/bin/env python3
import pysam
import random
from Bio.Seq import Seq
import pandas as pd
import subprocess
import os
from typing import List, Dict
import click
from tabulate import tabulate

def create_test_primers() -> pd.DataFrame:
    """Create test primer pairs"""
    primers = {
        'Name': ['Pair1', 'Pair2', 'Pair3'],
        'Forward': ['ACGTACGTACGT', 'GCTAGCTAGCTA', 'TGCATGCATGCA'],
        'Reverse': ['TGCATGCATGCA', 'TAGCTAGCTACG', 'ACGTACGTACGT'],
        'Size': [100, 150, 200]
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

def get_valid_size(size: int, primers_df: pd.DataFrame) -> int:
    """Get size within tolerance"""
    return size - len(primers_df.iloc[0]['Forward']) - len(primers_df.iloc[0]['Reverse'])

def get_invalid_size(size: int) -> int:
    """Get size outside tolerance"""
    return max(20, int(size * 0.5))  # Half the expected size, minimum 20bp

def generate_correct_pairs(i: int, row: pd.Series, primers_df: pd.DataFrame) -> str:
    return (row['Forward'] + 
            generate_read(get_valid_size(row['Size'], primers_df)) + 
            str(Seq(row['Reverse']).reverse_complement()))

def generate_wrong_size(i: int, row: pd.Series, primers_df: pd.DataFrame) -> str:
    return (row['Forward'] + 
            generate_read(get_invalid_size(row['Size'])) + 
            str(Seq(row['Reverse']).reverse_complement()))

def generate_single_end(i: int, row: pd.Series) -> str:
    return (row['Forward'] + 
            generate_read(get_invalid_size(row['Size'])))

def generate_mismatched_pairs(i: int, row: pd.Series, primers_df: pd.DataFrame) -> str:
    return (row['Forward'] + 
            generate_read(get_valid_size(row['Size'], primers_df)) + 
            str(Seq(primers_df.iloc[(i + 1) % len(primers_df)]['Reverse']).reverse_complement()))

def generate_one_primer_terminal(i: int, row: pd.Series, primers_df: pd.DataFrame, valid_size: bool = True) -> str:
    TERMINUS_LENGTH = 15
    size_func = get_valid_size if valid_size else get_invalid_size
    size = size_func(row['Size'], primers_df) if valid_size else size_func(row['Size'])
    return (row['Forward'] + 
            generate_read(size) + 
            str(Seq(row['Reverse']).reverse_complement())[:TERMINUS_LENGTH])

def generate_paired_terminal(i: int, row: pd.Series, primers_df: pd.DataFrame, valid_size: bool = True) -> str:
    TERMINUS_LENGTH = 15
    size_func = get_valid_size if valid_size else get_invalid_size
    size = size_func(row['Size'], primers_df) if valid_size else size_func(row['Size'])
    return (row['Forward'][:TERMINUS_LENGTH] + 
            generate_read(size) + 
            str(Seq(row['Reverse']).reverse_complement())[:TERMINUS_LENGTH])

def generate_single_terminal(i: int, row: pd.Series, primers_df: pd.DataFrame, valid_size: bool = True) -> str:
    TERMINUS_LENGTH = 15
    size_func = get_valid_size if valid_size else get_invalid_size
    size = size_func(row['Size'], primers_df) if valid_size else size_func(row['Size'])
    return (row['Forward'][:TERMINUS_LENGTH] + 
            generate_read(size))

def generate_multi_primer(i: int, row: pd.Series, primers_df: pd.DataFrame) -> str:
    return (row['Forward'] + row['Forward'] +  # Duplicate primer
            generate_read(get_valid_size(row['Size'], primers_df)) + 
            str(Seq(row['Reverse']).reverse_complement()))

def generate_reverse_orientation(i: int, row: pd.Series, primers_df: pd.DataFrame) -> str:
    return (str(Seq(row['Reverse']).reverse_complement()) + 
            generate_read(get_valid_size(row['Size'], primers_df)) + 
            row['Forward'])

def create_test_reads(categories: Dict = None) -> List[Dict]:
    """Create test reads for each category"""
    reads = []
    primers_df = create_test_primers()
    
    if categories is None:
        categories = {
            'no_primers': {
                'count': 400,  # ~35% of reads
                'category': '🟥 No primers or terminal matches detected',
                'generator': lambda i, row: generate_read(random.randint(40, 60))
            },
            'correct_pairs': {
                'count': 26,  # ~2% of reads
                'category': '🟩 Matched pairs - correct orientation and size',
                'generator': lambda i, row: generate_correct_pairs(i, row, primers_df)
            },
            'wrong_size': {
                'count': 523,  # ~47% of reads
                'category': '🟧 Matched pairs - correct orientation, wrong size',
                'generator': lambda i, row: generate_wrong_size(i, row, primers_df)
            },
            'single_end': {
                'count': 168,  # ~15% of reads
                'category': '🟨 Single-end primers only (no terminal match)',
                'generator': generate_single_end
            },
            'mismatched_pairs': {
                'count': 2,  # <1% of reads
                'category': '🟥 Mismatched primer pairs (different primers)',
                'generator': lambda i, row: generate_mismatched_pairs(i, row, primers_df)
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

@click.command()
@click.option('--output-dir', default='test_data', help='Output directory for test files')
@click.option('--verbose', is_flag=True, help='Show detailed output')
def main(output_dir: str, verbose: bool):
    """Generate test dataset and run URAdime analysis"""
    os.makedirs(output_dir, exist_ok=True)
    
    with click.progressbar(length=4, label='Generating test data') as bar:
        primers_df = create_test_primers()
        primers_df.to_csv(f'{output_dir}/test_primers.tsv', sep='\t', index=False)
        bar.update(1)
        
        reads = create_test_reads()
        bam_path = f'{output_dir}/test.bam'
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
        bar.update(1)
        
        # Run URAdime
        click.echo("Running URAdime analysis...")
        cmd = [
            'uradime',
            '-b', bam_path,
            '-p', f'{output_dir}/test_primers.tsv',
            '-o', f'{output_dir}/uradime_results',
            '--max-distance', '0'
        ]
        
        subprocess.run(cmd, capture_output=not verbose)
        bar.update(1)
        
        # Load actual results
        actual_df = pd.read_csv(f'{output_dir}/uradime_results_summary.csv')
        
        # Merge and display comparison
        comparison = pd.merge(
            expected_df, 
            actual_df, 
            on='Category', 
            suffixes=('_expected', '_actual'),
            how='outer'
        ).fillna(0)
        
        click.echo("\nTest Results:")
        click.echo(format_comparison_table(comparison))

if __name__ == '__main__':
    main()