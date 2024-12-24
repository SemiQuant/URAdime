#!/usr/bin/env python

import os
import sys
import pysam
from Bio import SeqIO

def fasta_to_bam(fasta_file, bam_file):
    """Convert FASTA file to BAM format"""
    # Create header
    header = {'HD': {'VN': '1.0'},
             'SQ': [{'LN': 1000, 'SN': 'test_ref'}]}
    
    # Open output BAM file
    with pysam.AlignmentFile(bam_file, 'wb', header=header) as outf:
        # Read sequences from FASTA
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Create alignment record
            a = pysam.AlignedSegment()
            a.query_name = record.id
            a.query_sequence = str(record.seq)
            a.flag = 4  # unmapped
            a.is_unmapped = True
            outf.write(a)
    
    # Sort and index BAM file
    pysam.sort('-o', bam_file + '.sorted', bam_file)
    os.rename(bam_file + '.sorted', bam_file)
    pysam.index(bam_file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.fasta output.bam")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    bam_file = sys.argv[2]
    fasta_to_bam(fasta_file, bam_file) 