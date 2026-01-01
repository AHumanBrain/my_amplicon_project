import csv
from Bio import SeqIO

def parse_bed(bed_path):
    """Parses a BED file into a list of (chrom, start, end, name) tuples."""
    intervals = []
    with open(bed_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                intervals.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return intervals

def parse_target_list(file_path):
    """Reads a list of target names from a text file."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def write_csv_recipe(output_path, primer_data):
    """Writes a standardized CSV for primer synthesis ordering."""
    headers = ['Sequence Name', 'Sequence', 'Target', 'Type', 'TM', 'Volume (uL)']
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(primer_data)
