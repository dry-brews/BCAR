#!/usr/bin/env python3
"""
FASTQ Barcode Splitter

Splits a sorted FASTQ file into smaller files based on unique barcode counts.
Assumes the input file is already sorted by barcode.
"""

import argparse
import sys
from pathlib import Path


def parse_fastq_record(lines):
    """Parse a single FASTQ record from 4 lines."""
    if len(lines) != 4:
        return None
    
    header = lines[0].strip()
    sequence = lines[1].strip()
    plus_line = lines[2].strip()
    quality = lines[3].strip()
    
    return {
        'header': header,
        'sequence': sequence,
        'plus_line': plus_line,
        'quality': quality
    }


def write_fastq_record(file_handle, record):
    """Write a FASTQ record to file."""
    file_handle.write(f"{record['header']}\n")
    file_handle.write(f"{record['sequence']}\n")
    file_handle.write(f"{record['plus_line']}\n")
    file_handle.write(f"{record['quality']}\n")


def extract_barcode(sequence, start_pos, length):
    """Extract barcode from sequence at specified position and length."""
    end_pos = start_pos + length
    if end_pos > len(sequence):
        return None
    return sequence[start_pos:end_pos]


def generate_output_filename(input_path, file_number):
    """Generate output filename with incrementing number."""
    input_path = Path(input_path)
    stem = input_path.stem
    suffix = input_path.suffix
    
    # Handle compressed files (.fastq.gz)
    if suffix == '.gz' and stem.endswith('.fastq'):
        stem = stem[:-6]  # Remove .fastq
        suffix = '.fastq.gz'
    elif suffix == '.gz' and stem.endswith('.fq'):
        stem = stem[:-3]  # Remove .fq
        suffix = '.fq.gz'
    
    output_filename = f"{stem}_{file_number:04d}{suffix}"
    return input_path.parent / output_filename


def split_fastq_by_barcodes(input_file, barcode_start, barcode_length, 
                           barcodes_per_file, output_dir=None):
    """
    Split FASTQ file based on unique barcode counts.
    
    Args:
        input_file: Path to input FASTQ file
        barcode_start: Starting position of barcode in sequence (0-based)
        barcode_length: Length of barcode
        barcodes_per_file: Number of unique barcodes per output file
        output_dir: Output directory (default: same as input file)
    """
    
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    output_dir = input_path.parent
    
    # Initialize variables
    current_file_num = 1
    current_barcode_count = 0
    current_barcode = None
    current_output_file = None
    total_reads = 0
    total_barcodes = 0
    
    print(f"Processing: {input_file}")
    print(f"Barcode position: {barcode_start}, length: {barcode_length}")
    print(f"Barcodes per file: {barcodes_per_file}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    try:
        # Handle both regular and gzipped files
        if input_file.endswith('.gz'):
            import gzip
            file_opener = gzip.open
            mode = 'rt'
        else:
            file_opener = open
            mode = 'r'
        
        with file_opener(input_file, mode) as infile:
            lines = []
            
            for line in infile:
                lines.append(line)
                
                # Process every 4 lines (one FASTQ record)
                if len(lines) == 4:
                    record = parse_fastq_record(lines)
                    if record is None:
                        print(f"Warning: Invalid FASTQ record at read {total_reads + 1}")
                        lines = []
                        continue
                    
                    # Extract barcode
                    barcode = extract_barcode(record['sequence'], barcode_start, barcode_length)
                    if barcode is None:
                        print(f"Warning: Could not extract barcode from read {total_reads + 1}")
                        lines = []
                        continue
                    
                    # Check if we need to start a new output file
                    if current_output_file is None:
                        # First file
                        output_filename = generate_output_filename(input_path, current_file_num)
                        output_path = output_dir / output_filename.name
                        
                        if output_path.suffix.endswith('.gz'):
                            current_output_file = gzip.open(output_path, 'wt')
                        else:
                            current_output_file = open(output_path, 'w')
                        
                        current_barcode = barcode
                        current_barcode_count = 1
                        total_barcodes = 1
                        print(f"Started file {current_file_num}: {output_path.name}")
                    
                    elif barcode != current_barcode:
                        # New barcode encountered
                        current_barcode_count += 1
                        current_barcode = barcode
                        total_barcodes += 1
                        
                        # Check if we need to start a new file
                        if current_barcode_count > barcodes_per_file:
                            current_output_file.close()
                            current_file_num += 1
                            current_barcode_count = 1
                            
                            output_filename = generate_output_filename(input_path, current_file_num)
                            output_path = output_dir / output_filename.name
                            
                            if output_path.suffix.endswith('.gz'):
                                current_output_file = gzip.open(output_path, 'wt')
                            else:
                                current_output_file = open(output_path, 'w')
                            
                            print(f"Started file {current_file_num}: {output_path.name}")
                    
                    # Write record to current output file
                    write_fastq_record(current_output_file, record)
                    total_reads += 1
                    
                    # Progress reporting
                    if total_reads % 100000 == 0:
                        print(f"Processed {total_reads:,} reads, {total_barcodes:,} unique barcodes")
                    
                    lines = []
        
        # Close the last output file
        if current_output_file:
            current_output_file.close()
        
        print("-" * 50)
        print("Processing complete!")
        print(f"Total reads processed: {total_reads:,}")
        print(f"Total unique barcodes: {total_barcodes:,}")
        print(f"Output files created: {current_file_num}")
        
    except Exception as e:
        if current_output_file:
            current_output_file.close()
        print(f"Error: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Split sorted FASTQ file based on unique barcode counts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split file with barcodes at position 0, length 8, 1000 barcodes per file
  python fastq_splitter.py input.fastq --start 0 --length 8 --count 1000
  
  # Split with custom output directory
  python fastq_splitter.py input.fastq.gz --start 10 --length 12 --count 500 --output-dir split_files/
        """
    )
    
    parser.add_argument('--in1', 
                       type=str, required=True,
                       help='Input FASTQ file (can be gzipped)')
    
    parser.add_argument('--bc-start', '-s', 
                       type=int, required=True,
                       help='Starting position of barcode in sequence (0-based)')
    
    parser.add_argument('--bc-len', '-l', 
                       type=int, required=True,
                       help='Length of barcode')
    
    parser.add_argument('--size', '-n', 
                       type=int, required=True,
                       help='Number of unique barcodes per output file')

    args = parser.parse_args()
    
    # Validate arguments
    if args.bc_start < 0:
        print("Error: Barcode start position must be non-negative")
        sys.exit(1)
    
    if args.bc_len <= 0:
        print("Error: Barcode length must be positive")
        sys.exit(1)
    
    if args.size <= 0:
        print("Error: Barcode count per file must be positive")
        sys.exit(1)
    
    # Run the splitter
    split_fastq_by_barcodes(
        args.in1,
        args.bc_start,
        args.bc_len,
        args.size,
    )


if __name__ == "__main__":
    main()