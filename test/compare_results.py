import click
import pandas as pd
import pysam
import re
from tabulate import tabulate
from collections import defaultdict
from pathlib import Path

def analyze_toy_dataset(bam_file):
    """Analyze the toy dataset BAM file to get ground truth statistics."""
    read_patterns = {
        'perfect_match': r'^perfect_match_(\w+)_\d+$',
        'wrong_size': r'^wrong_size_(\w+)_\d+$',
        'single_terminal': r'^single_terminal_(start|end)_(\w+)_\d+$',
        'paired_terminal': r'^paired_terminal_(\w+)_\d+$',
        'single_end': r'^single_end_(forward|reverse)_(\w+)_\d+$',
        'mismatched_pair': r'^mismatched_pair_(\w+)_\d+$',
        'no_primers': r'^no_primers_(\w+)_\d+$',
        'hybrid': r'^hybrid_(forward|reverse)_full_(\w+)_\d+$'
    }
    
    stats = defaultdict(lambda: defaultdict(int))
    total_reads = 0
    
    with click.progressbar(pysam.AlignmentFile(bam_file, 'rb'),
                         label='Analyzing toy dataset') as bam:
        for read in bam:
            total_reads += 1
            read_name = read.query_name
            
            for pattern_type, pattern in read_patterns.items():
                match = re.match(pattern, read_name)
                if match:
                    if pattern_type in ['single_terminal', 'hybrid']:
                        primer_name = match.group(2)
                    elif pattern_type == 'single_end':
                        primer_name = match.group(2)
                    else:
                        primer_name = match.group(1)
                    
                    stats[primer_name][pattern_type] += 1
                    break
    
    return stats, total_reads

def analyze_uradime_output(summary_file, matched_file, mismatched_file, wrong_size_file):
    """Analyze URADime output files."""
    stats = defaultdict(lambda: defaultdict(int))
    categories = {}
    
    with click.progressbar(length=4, label='Analyzing URADime output') as bar:
        # Process matched pairs (perfect matches)
        if Path(matched_file).exists():
            matched_df = pd.read_csv(matched_file)
            for _, row in matched_df.iterrows():
                primer_name = row['Start_Primer_Name']
                if row['Correct_Orientation'] and row['Size_Compliant']:
                    stats[primer_name]['perfect_matches'] += 1
        bar.update(1)
        
        # Process wrong size pairs
        if Path(wrong_size_file).exists():
            wrong_size_df = pd.read_csv(wrong_size_file)
            for _, row in wrong_size_df.iterrows():
                primer_name = row['Start_Primer_Name']
                stats[primer_name]['wrong_size'] += 1
        bar.update(1)
        
        # Process mismatched pairs
        if Path(mismatched_file).exists():
            mismatched_df = pd.read_csv(mismatched_file)
            for _, row in mismatched_df.iterrows():
                primer_name = row['Start_Primer_Name']
                stats[primer_name]['mismatched'] += 1
        bar.update(1)
        
        # Process summary
        if Path(summary_file).exists():
            summary_df = pd.read_csv(summary_file)
            # Extract categories from summary
            for _, row in summary_df.iterrows():
                categories[row['Category']] = {
                    'Count': row['Count'],
                    'Percentage': row['Percentage']
                }
            
            # Group by primer name and calculate statistics
            for primer_name in stats.keys():
                total_reads = 0
                correct_orientation = 0
                size_compliant = 0
                
                # Sum up statistics from all files
                if Path(matched_file).exists():
                    primer_reads = matched_df[matched_df['Start_Primer_Name'] == primer_name]
                    total_reads += len(primer_reads)
                    correct_orientation += primer_reads['Correct_Orientation'].sum()
                    size_compliant += primer_reads['Size_Compliant'].sum()
                
                if Path(wrong_size_file).exists():
                    primer_reads = wrong_size_df[wrong_size_df['Start_Primer_Name'] == primer_name]
                    total_reads += len(primer_reads)
                    correct_orientation += primer_reads['Correct_Orientation'].sum()
                
                if Path(mismatched_file).exists():
                    primer_reads = mismatched_df[mismatched_df['Start_Primer_Name'] == primer_name]
                    total_reads += len(primer_reads)
                
                if total_reads > 0:
                    stats[primer_name]['total'] = total_reads
                    stats[primer_name]['correct_orientation'] = (correct_orientation / total_reads) * 100
                    stats[primer_name]['size_compliant'] = (size_compliant / total_reads) * 100
        bar.update(1)
    
    return stats, categories

def compare_results(toy_stats, uradime_stats):
    """Compare toy dataset ground truth with URADime results."""
    # First table: Per-primer comparison
    primer_comparison = []
    stats, categories = uradime_stats
    
    for primer_name in sorted(set(toy_stats.keys()) | set(stats.keys())):
        toy = toy_stats[primer_name]
        ura = stats[primer_name]
        
        row = {
            'Primer': primer_name,
            'Perfect (Toy/URA)': f"{toy['perfect_match']}/{ura.get('perfect_matches', 0)}",
            'Wrong Size (Toy/URA)': f"{toy['wrong_size']}/{ura.get('wrong_size', 0)}",
            'Single End': toy['single_end'],
            'Mismatched (Toy/URA)': f"{toy['mismatched_pair']}/{ura.get('mismatched', 0)}",
            'Hybrid': toy['hybrid'],
            'No Primers': toy['no_primers'],
            'Total (Toy/URA)': f"{sum(toy.values())}/{ura.get('total', 0)}",
            'Orient.%': f"{ura.get('correct_orientation', 0):.0f}%",
            'Size%': f"{ura.get('size_compliant', 0):.0f}%"
        }
        primer_comparison.append(row)
    
    # Create a mapping between toy dataset and URAdime categories
    category_mapping = {
        'perfect_match': '🟩 Matched pairs - correct orientation and size',
        'wrong_size': '🟧 Matched pairs - correct orientation, wrong size',
        'single_terminal': '🟨 Single terminal match only - correct size',
        'paired_terminal': '🟨 Paired terminal matches - correct size',
        'single_end': '🟨 Single-end primers only (no terminal match)',
        'mismatched_pair': '🟥 Mismatched primer pairs (different primers)',
        'no_primers': '🟥 No primers or terminal matches detected',
        'hybrid': '🟨 One full primer + one terminal match - correct size'
    }
    
    # Combine toy dataset and URAdime results
    comparison_rows = []
    total_toy_reads = sum(sum(primer_stats.values()) for primer_stats in toy_stats.values())
    matches = 0
    total_categories = 0
    
    # Add categories from toy dataset
    for toy_category, uradime_category in category_mapping.items():
        toy_count = sum(primer_stats[toy_category] for primer_stats in toy_stats.values())
        uradime_data = categories.get(uradime_category, {'Count': 0, 'Percentage': 0.0})
        is_match = toy_count == uradime_data['Count']
        
        if toy_count > 0 or uradime_data['Count'] > 0:
            total_categories += 1
            if is_match:
                matches += 1
        
        row = {
            'Category': uradime_category,
            'Expected Count': toy_count,
            'Expected %': f"{(toy_count / total_toy_reads * 100):.1f}%",
            'URAdime Count': uradime_data['Count'],
            'URAdime %': f"{uradime_data['Percentage']:.1f}%",
            'Match': '✅' if is_match else '❌'
        }
        comparison_rows.append(row)
    
    # Add any remaining URAdime categories not in toy dataset
    for category, data in categories.items():
        if category not in category_mapping.values():
            if data['Count'] > 0:
                total_categories += 1
            row = {
                'Category': category,
                'Expected Count': 0,
                'Expected %': '0.0%',
                'URAdime Count': data['Count'],
                'URAdime %': f"{data['Percentage']:.1f}%",
                'Match': '❌'
            }
            comparison_rows.append(row)
    
    category_df = pd.DataFrame(comparison_rows)
    primer_df = pd.DataFrame(primer_comparison)
    
    # Sort categories to match URAdime's order
    category_order = {cat: i for i, cat in enumerate(categories.keys())}
    category_df['sort_order'] = category_df['Category'].map(category_order)
    category_df = category_df.sort_values('sort_order').drop('sort_order', axis=1)
    
    # Calculate overall agreement percentage
    agreement_percentage = (matches / total_categories * 100) if total_categories > 0 else 0
    
    return primer_df, category_df, agreement_percentage 