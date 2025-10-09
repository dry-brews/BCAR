#!/usr/bin/env python3
"""
Compare two FASTQ files to determine if they contain identical reads.
Reads may be in different orders.
"""

import sys
import argparse
from collections import defaultdict


def parse_fastq(filename):
    """
    Parse a FASTQ file and return a set of read tuples.
    Each tuple contains (header, sequence, plus_line, quality).
    """
    reads = set()
    
    with open(filename, 'r') as f:
        while True:
            header = f.readline().rstrip()
            if not header:
                break
            
            sequence = f.readline().rstrip()
            plus_line = f.readline().rstrip()
            quality = f.readline().rstrip()
            
            # Store as tuple (immutable and hashable)
            reads.add((header, sequence, plus_line, quality))
    
    return reads


def compare_fastq_files(file1, file2):
    """
    Compare two FASTQ files for identical reads.
    Returns True if files contain identical reads, False otherwise.
    """
    try:
        reads1 = parse_fastq(file1)
        reads2 = parse_fastq(file2)
        
        return reads1 == reads2
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Compare two FASTQ files for identical reads'
    )
    parser.add_argument('file1', help='First FASTQ file')
    parser.add_argument('file2', help='Second FASTQ file')
    
    args = parser.parse_args()
    
    result = compare_fastq_files(args.file1, args.file2)
    print(result)


if __name__ == '__main__':
    main()
