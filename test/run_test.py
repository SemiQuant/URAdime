import click
import subprocess
import os
import sys
from pathlib import Path
from tabulate import tabulate

# Add the test directory to the Python path
test_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(test_dir)

from compare_results import analyze_toy_dataset, analyze_uradime_output, compare_results

def format_table(df, max_width=120):
    """Format DataFrame into a nice table with wrapped text."""
    # Rename columns to be more concise
    column_renames = {
        'Expected Count': 'Exp.Count',
        'Expected %': 'Exp.%',
        'URAdime Count': 'URA.Count',
        'URAdime %': 'URA.%'
    }
    df = df.rename(columns=column_renames)
    
    # Format the table
    table = tabulate(
        df,
        headers='keys',
        tablefmt='pretty',
        showindex=False,
        numalign='right',
        stralign='left'
    )
    
    # Split table into lines and ensure no line is longer than max_width
    lines = table.split('\n')
    formatted_lines = []
    for line in lines:
        if len(line) > max_width:
            # Wrap the line
            current_pos = 0
            wrapped_line = ''
            while current_pos < len(line):
                wrapped_line += line[current_pos:current_pos + max_width] + '\n'
                current_pos += max_width
            formatted_lines.append(wrapped_line.rstrip())
        else:
            formatted_lines.append(line)
    
    return '\n'.join(formatted_lines)

def write_report(output_dir, category_df, primer_df, agreement_percentage):
    """Write a detailed report to a file."""
    report_path = os.path.join(output_dir, 'test_report.txt')
    with open(report_path, 'w') as f:
        f.write("URAdime Test Report\n")
        f.write("==================\n\n")
        
        f.write("Overall Agreement: {:.1f}%\n\n".format(agreement_percentage))
        
        f.write("Category Comparison:\n")
        f.write("-------------------\n")
        f.write(format_table(category_df))
        f.write("\n\n")
        
        f.write("Per-Primer Comparison:\n")
        f.write("---------------------\n")
        f.write(format_table(primer_df))
        f.write("\n")

@click.command()
@click.option('--output-dir', default='test/results', help='Directory to store URAdime results')
def main(output_dir):
    """Run URAdime test workflow:
    1. Generate toy dataset
    2. Run URAdime with max-distance 0
    3. Compare results with expected data
    """
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    with click.progressbar(length=3, label='Running test workflow') as bar:
        # Step 1: Generate toy dataset
        click.echo("\n🧬 Generating toy dataset...")
        subprocess.run(['python', 'test/generate_toy_dataset.py'], check=True)
        bar.update(1)
        
        # Step 2: Run URAdime
        click.echo("\n🔍 Running URAdime analysis...")
        uradime_prefix = os.path.join(output_dir, 'uradime')
        subprocess.run([
            'uradime',
            '--max-distance', '0',
            '--bam', 'test/data/test_data.bam',
            '--primers', 'test/data/test_primers.tsv',
            '-o', uradime_prefix
        ], check=True)
        bar.update(1)
        
        # Step 3: Compare results
        click.echo("\n📊 Comparing results...")
        toy_stats, total_toy_reads = analyze_toy_dataset('test/data/test_data.bam')
        uradime_stats = analyze_uradime_output(
            f"{uradime_prefix}_summary.csv",
            f"{uradime_prefix}_matched_pairs.csv",
            f"{uradime_prefix}_mismatched_pairs.csv",
            f"{uradime_prefix}_wrong_size_pairs.csv"
        )
        
        primer_df, category_df, agreement_percentage = compare_results(toy_stats, uradime_stats)
        
        # Save comparison results
        primer_df.to_csv(os.path.join(output_dir, 'primer_comparison.csv'), index=False)
        category_df.to_csv(os.path.join(output_dir, 'category_summary.csv'), index=False)
        write_report(output_dir, category_df, primer_df, agreement_percentage)
        
        # Display results in nice tables
        click.echo("\n📊 Overall Agreement: {:.1f}%".format(agreement_percentage))
        
        click.echo("\n📋 URAdime Categories Summary:")
        click.echo(format_table(category_df))
        
        # click.echo("\n📋 Per-Primer Comparison:")
        # click.echo(format_table(primer_df))
        # bar.update(1)

if __name__ == "__main__":
    main() 