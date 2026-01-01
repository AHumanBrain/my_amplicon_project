import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import os

def measure_uniqueness(fasta_path, start_k=12, end_k=25):
    """
    Measures the percentage of non-unique k-mers for a range of k values.
    """
    print(f"Analyzing genome: {fasta_path}")
    
    # Load Genome
    genome_seq = ""
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seq += str(record.seq).upper()
    
    print(f"Genome Length: {len(genome_seq)} bp")
    
    results_k = []
    results_pct = []

    print(f"{'K-mer':<10} | {'Total':<15} | {'Non-Unique':<15} | {'% Non-Unique':<15}")
    print("-" * 65)

    for k in range(start_k, end_k + 1):
        kmer_counts = {}
        total_kmers = 0
        
        # Sliding window
        for i in range(len(genome_seq) - k + 1):
            kmer = genome_seq[i:i+k]
            # Simple optimization: only count canonical k-mers? 
            # ideally yes, but for 'collision' check, direct string exact match is sufficient proxy
            if 'N' in kmer: continue
            
            total_kmers += 1
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
        
        # Calculate statistics
        non_unique_count = sum(1 for count in kmer_counts.values() if count > 1)
        # Wait, "percent of primers that are non-unique". 
        # Does this mean "Percent of unique sequences that appear >1 time"?
        # Or "Percent of the genome covered by non-unique kmers"?
        # User asked: "percent of primers... that are non-unique". 
        # Interpretation: If I pick a random k-mer from the genome, what is the chance it appears elsewhere?
        # This corresponds to: sum(counts where count > 1) / total_kmers
        
        non_unique_instances = sum(count for count in kmer_counts.values() if count > 1)
        pct_non_unique = (non_unique_instances / total_kmers) * 100.0
        
        results_k.append(k)
        results_pct.append(pct_non_unique)
        
        print(f"{k:<10} | {total_kmers:<15} | {non_unique_instances:<15} | {pct_non_unique:.2f}%")

    return results_k, results_pct

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python measure_kmer_uniqueness.py <path_to_fasta>")
        sys.exit(1)
        
    fasta_file = sys.argv[1]
    ks, pcts = measure_uniqueness(fasta_file)
    
    # Plotting (headless)
    try:
        plt.figure(figsize=(10, 6))
        plt.plot(ks, pcts, marker='o', linestyle='-', color='b')
        plt.title(f"Genomic Non-Uniqueness vs K-mer Length\nFile: {os.path.basename(fasta_file)}")
        plt.xlabel("K-mer Length (bp)")
        plt.ylabel("% Non-Unique Hit Probability")
        plt.grid(True)
        plt.xticks(ks)
        plt.ylim(0, 100)
        
        output_plot = "kmer_uniqueness_plot.png"
        plt.savefig(output_plot)
        print(f"\nPlot saved to {output_plot}")
    except Exception as e:
        print(f"Could not generate plot: {e}")
