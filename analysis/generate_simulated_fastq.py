import argparse
import gzip
import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_genome(genome_path):
    print(f"Loading genome from {genome_path}...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    return genome_dict

def parse_bed(bed_path):
    print(f"Parsing BED file {bed_path}...")
    intervals = []
    with open(bed_path, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                intervals.append((chrom, start, end, name))
    return intervals

def mutate_sequence(sequence, mutation_rate=0.01):
    """Introduces random SNVs into the sequence."""
    seq_list = list(sequence)
    bases = ['A', 'C', 'G', 'T']
    mutations = []
    
    for i in range(len(seq_list)):
        if random.random() < mutation_rate:
            original_base = seq_list[i]
            possible_bases = [b for b in bases if b != original_base]
            new_base = random.choice(possible_bases)
            seq_list[i] = new_base
            mutations.append((i, original_base, new_base))
            
    return "".join(seq_list), mutations

def generate_fastq_records(intervals, genome_dict, coverage, read_len, output_prefix, non_uniform=False, thermal_bias=None):
    r1_records = []
    r2_records = []
    
    total_reads = 0
    
    for chrom, start, end, name in intervals:
        if chrom not in genome_dict:
            print(f"Warning: Chromosome {chrom} not found in genome. Skipping {name}.")
            continue
            
        ref_seq = str(genome_dict[chrom].seq[start:end])
        
        # Inject mutations to simulate biological variation
        mutated_seq, mutations = mutate_sequence(ref_seq, mutation_rate=0.01)
        if mutations:
           print(f"   -> Injected {len(mutations)} SNVs into {name}")

        # Determine num_pairs for this amplicon
        current_coverage = coverage
        
        if thermal_bias:
            # Simple simulation: amplicons have an 'optimal' temp around 62C
            # Coverage drops off as thermal_bias moves away from optimal
            # We'll use the amplicon name or a hash to vary the 'optimal' temp slightly per-target
            target_opt = 62.0 + (hash(name) % 5) - 2.5 # 60 to 65C range
            diff = abs(thermal_bias - target_opt)
            # Scaling factor: 1.0 at opt, drops by 0.2 for every degree of diff
            scaling = max(0.05, 1.0 - (diff * 0.2))
            current_coverage = int(coverage * scaling)

        if non_uniform:
            # Randomly vary coverage between 10% and 200% of target
            num_pairs = int(current_coverage * random.uniform(0.1, 2.0))
        else:
            num_pairs = current_coverage
        
        print(f"   -> Generating {num_pairs} pairs for {name} (Bias: {thermal_bias}C)")
        
        for i in range(num_pairs):
            # R1: Forward strand, 5' end
            r1_seq = mutated_seq[:read_len]
            if len(r1_seq) < read_len:
                r1_seq = r1_seq.ljust(read_len, 'N')
                
            # R2: Reverse strand, 3' end (reverse complement)
            r2_seq_template = mutated_seq[-read_len:]
            if len(r2_seq_template) < read_len:
                r2_seq_template = r2_seq_template.rjust(read_len, 'N')
            
            r2_seq = str(Seq(r2_seq_template).reverse_complement())
            qual = "I" * read_len
            header = f"SimonSimRead:{total_reads}:{name}"
            
            r1_records.append(f"@{header} 1:N:0:1\n{r1_seq}\n+\n{qual}\n")
            r2_records.append(f"@{header} 2:N:0:1\n{r2_seq}\n+\n{qual}\n")
            total_reads += 1

    r1_filename = f"{output_prefix}_R1.fastq.gz"
    r2_filename = f"{output_prefix}_R2.fastq.gz"

    print(f"Writing {len(r1_records)} read pairs to {r1_filename} and {r2_filename}...")

    with gzip.open(r1_filename, 'wt') as f1, gzip.open(r2_filename, 'wt') as f2:
        for r1, r2 in zip(r1_records, r2_records):
            f1.write(r1)
            f2.write(r2)

def main():
    parser = argparse.ArgumentParser(description="Generate simulated FASTQ data from BED intervals.")
    parser.add_argument('--genome', required=True, help="Path to reference genome")
    parser.add_argument('--bed', required=True, help="Path to BED file")
    parser.add_argument('--output-dir', required=True, help="Directory to save FASTQ files")
    parser.add_argument('--sample-name', default="Simulated_Sample", help="Sample name for output files")
    parser.add_argument('--coverage', type=int, default=100, help="Number of read pairs per amplicon")
    parser.add_argument('--non-uniform', action='store_true', help="Introduce coverage non-uniformity")
    parser.add_argument('--thermal-bias', type=float, help="Simulated annealing temperature (C). Affects coverage of amplicons based on proximity to 60C.")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    output_prefix = os.path.join(args.output_dir, args.sample_name)
    
    genome = load_genome(args.genome)
    intervals = parse_bed(args.bed)
    
    generate_fastq_records(intervals, genome, args.coverage, 150, output_prefix, 
                           non_uniform=args.non_uniform, thermal_bias=args.thermal_bias)
    print("Simulation complete.")

if __name__ == "__main__":
    main()
