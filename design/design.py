import argparse
from Bio import SeqIO
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq # Added for reverse complement
import primer3
import os
import csv
import subprocess
import shutil
import multiprocessing
import itertools
from functools import partial, lru_cache
from tqdm import tqdm
import re # Import re for parsing

# --- (v7.7) Refactored Hardcoded Adapter Tails ---
# Tails are the full length from index to insert
#FWD_P5_TRUNCATED = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT' 33-mer
#REV_P7_TRUNCATED = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' 34-mer

# --- (v8.0 Ready) Refactored Hardcoded Adapter Tails ---
# Tails shortened to 20bp to keep final oligo length < 60nt
# This allows for standard desalting without loss of the 5' clamp.

# For 'fwd_tailed' (synthesis-ready) format
FWD_P5_TRUNCATED = 'ACACGACGCTCTTCCGATCT' #20-mer; shaved off the 5'-most 13 bases
REV_P7_TRUNCATED = 'GACGTGTGCTCTTCCGATCT' #20-mer; shaved off the 5' most 14 bases

# For 'rc_tailed' (template-generation) format
# These remain long as they are not for direct synthesis.
FWD_RC_P5_TRUNCATED = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
REV_RC_P7_TRUNCATED = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# --- Design Constants ---
END_STABILITY_DG_THRESHOLD = -9.0
IDEAL_TM_MIN = 59.0
IDEAL_TM_MAX = 61.0
# (v8.0 Ready) Created a single source of truth for min primer size
IDEAL_PRIMER_MIN_SIZE = 19 

# --- (v8.0 Ready) Phase 2 Thermodynamic Thresholds ---
# Lower (more negative) means stronger/worse dimerization
MAX_HETERODIMER_DG = -8.0  # (kcal/mol) Stop heterodimers at internal secondary structures
MAX_SELF_DIMER_DG = -9.0   # (kcal/mol) Stop primers folding back on themselves

# --- Large-Panel Logic ---
MAX_CLASH_RECOMMENDATION = 5
# (v8.0 Ready) Set the default for iterations to 5000
MAX_COMPATIBILITY_ITERATIONS = 5000

# --- (v8.0 Ready) Constants for Hairpin Clamp Logic ---
# Calibrated based on v7.5.1 results.
HAIRPIN_STEM_TARGET_DG = -11.9 # (kcal/mol) Calibrated target dG for the *stem* interaction
HAIRPIN_CLAMP_MAX_LEN = IDEAL_PRIMER_MIN_SIZE  # (v8.0 Ready) Linked to min primer size

# --- Core Helper Functions ---

def read_lines_from_file(file_path):
    """Reads a list of items from a file, one per line, skipping empty lines."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# --- (v7.8) Corrected Helper Function for Hairpin Clamps ---

# --- (v8.0 Ready) Phase 2 Cached Thermodynamic Functions ---
@lru_cache(maxsize=10000)
def get_heterodimer_dg(seq1, seq2):
    """Calculates worst-case dG across the entire sequence length."""
    try:
        return primer3.calc_heterodimer(seq1, seq2).dg / 1000.0
    except:
        return 0.0

@lru_cache(maxsize=10000)
def get_end_stability_dg(seq1, seq2):
    """Calculates stability of the 3' end interaction."""
    try:
        return primer3.calc_end_stability(seq1, seq2).dg / 1000.0
    except:
        return 0.0

@lru_cache(maxsize=10000)
def get_homodimer_dg(seq):
    """Calculates stability of a primer folding back on itself."""
    try:
        return primer3.calc_homodimer(seq).dg / 1000.0
    except:
        return 0.0

def add_iterative_hairpin_clamp(oligo_sequence, primer_specific_seq, target_stem_dg, max_clamp_len):
    """
    (v7.8) Iteratively adds clamp bases to the 5' end.
    Uses calc_heterodimer to model the stem interaction, aiming for the
    calibrated HAIRPIN_STEM_TARGET_DG.
    """
    try:
        primer_specific_seq_rc = str(Seq(primer_specific_seq).reverse_complement())
        
        best_oligo = oligo_sequence
        best_stem_dg = 0.0 # dG for 0-length clamp

        for i in range(1, min(max_clamp_len, len(primer_specific_seq_rc)) + 1):
            
            # (v7.8) Use [:i] to get complement of the 3' end.
            clamp_seq = primer_specific_seq_rc[:i]
            
            stem_dg = primer3.calc_heterodimer(clamp_seq, primer_specific_seq).dg / 1000.0

            if abs(stem_dg - target_stem_dg) < abs(best_stem_dg - target_stem_dg):
                best_oligo = clamp_seq + oligo_sequence
                best_stem_dg = stem_dg
            elif stem_dg < best_stem_dg: 
                break
            
        return best_oligo
    
    except Exception as e:
        print(f"Warning: Could not add hairpin clamp to {oligo_sequence[:20]}... Error: {e}")
        return oligo_sequence 

# --- Mode 1: Design Pipeline Functions ---

def create_blast_db_if_needed(genome_fasta_path, blast_db_prefix):
    db_files_exist = os.path.exists(f"{blast_db_prefix}.nin") or os.path.exists(f"{blast_db_prefix}.nhr")
    
    # Determine the cleaned filename
    base, ext = os.path.splitext(genome_fasta_path)
    if base.endswith('.cleaned'):
        cleaned_fasta_path = genome_fasta_path
    else:
        cleaned_fasta_path = f"{base}.cleaned{ext}"

    # Verify if we need to clean the genome
    cleaning_needed = not os.path.exists(cleaned_fasta_path)
    
    if db_files_exist and not cleaning_needed:
        print(f"BLAST database '{blast_db_prefix}' and cleaned genome already exist. Skipping creation.")
        return

    if not shutil.which("makeblastdb"):
        raise RuntimeError("`makeblastdb` command not found. Please ensure NCBI BLAST+ is installed and in your system's PATH.")

    if cleaning_needed:
        print(f"Cleaned genome file not found. Creating '{cleaned_fasta_path}' from '{genome_fasta_path}'...")
        print("   (This standardizes headers for BLAST and Primer3 compatibility)")
        
        cleaned_records = []
        try:
            for record in SeqIO.parse(genome_fasta_path, "fasta"):
                # Simplify header to just the ID to avoid BLAST/GFF mismatch issues
                record.id = record.id.split()[0]
                record.description = '' 
                cleaned_records.append(record)
        except Exception as e:
             raise ValueError(f"Biopython could not parse '{genome_fasta_path}'. It may not be a valid FASTA file. Original error: {e}")

        if not cleaned_records:
            raise ValueError(f"Input FASTA file '{genome_fasta_path}' is empty or contains no valid records. Please check the file format.")

        with open(cleaned_fasta_path, "w", newline='\n') as out_handle:
            SeqIO.write(cleaned_records, out_handle, "fasta")
    else:
        print(f"Using existing cleaned genome: '{cleaned_fasta_path}'")
    
    if not db_files_exist:
        print(f"Building BLAST database '{blast_db_prefix}'...")
        command = ['makeblastdb', '-in', cleaned_fasta_path, '-dbtype', 'nucl', '-out', blast_db_prefix]
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            print("BLAST database created successfully.")
        except subprocess.CalledProcessError as e:
            print("--- MAKEBLASTDB FAILED ---")
            print(f"Error: {e.stderr}")
            raise

def parse_gff(gff_file):
    gene_coords = {}
    gene_attr_re = re.compile(r"(?:Name|gene|locus_tag)=([^;]+)")

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene': continue
            
            contig, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attributes = parts[8]
            
            matches = gene_attr_re.findall(attributes)
            if not matches:
                continue

            for gene_name in matches:
                if gene_name not in gene_coords: 
                    gene_coords[gene_name] = {
                        'contig': contig, 
                        'start': start, 
                        'end': end, 
                        'strand': strand
                    }
                    
    if not gene_coords:
        print("Warning: GFF parsing finished but found 0 gene entries. Check GFF format and attributes (e.g., 'Name=', 'gene=', 'locus_tag=').")
        
    return gene_coords

def extract_target_sequence(genome_records, gene_coords, target_id, offset_bp=0):
    if target_id not in gene_coords:
        return None, None, f"Warning: Target ID '{target_id}' not found in GFF file. Skipping."
    
    target_info = gene_coords[target_id]
    contig_id = target_info['contig']
    
    if contig_id not in genome_records:
        core_contig_id = contig_id.split('|')[-1] if '|' in contig_id else contig_id
        if core_contig_id not in genome_records:
             return None, None, f"Warning: Contig '{contig_id}' (or '{core_contig_id}') for gene '{target_id}' not found in FASTA file. Skipping."
        target_info['contig'] = core_contig_id
    
    contig_seq = genome_records[target_info['contig']].seq
    start = max(0, target_info['start'] - 1 + offset_bp)
    end = min(len(contig_seq), target_info['end'] + offset_bp)
    target_seq = contig_seq[start:end]
    return str(target_seq), target_info, None

def design_primers_for_sequence(sequence, target_id, strategy_settings, num_candidates=25):
    seq_args = {'SEQUENCE_ID': target_id, 'SEQUENCE_TEMPLATE': sequence}
    
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': IDEAL_PRIMER_MIN_SIZE, # (v8.0 Ready) Use constant
        'PRIMER_MAX_SIZE': 21,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': IDEAL_TM_MIN,
        'PRIMER_MAX_TM': IDEAL_TM_MAX,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]],
        'PRIMER_NUM_RETURN': num_candidates
    }
    
    global_args.update(strategy_settings)
    
    return primer3.design_primers(seq_args, global_args)

def run_blast_specificity_check(primer_seq, blast_db_path, strict_threshold=1):
    """
    Check primer specificity using BLAST.
    
    Returns:
        tuple: (hit_count, hit_locations, is_specific, warning_message)
        - hit_count: Number of perfect matches in the genome
        - hit_locations: List of (contig, start, end) tuples
        - is_specific: True if hit_count <= strict_threshold
        - warning_message: Warning string if primer is non-specific, else None
    """
    command = [
        'blastn', '-query', '-', '-db', blast_db_path, '-task', 'blastn-short',
        '-outfmt', '6 sseqid sstart send', '-perc_identity', '100', '-qcov_hsp_perc', '100'
    ]
    try:
        process = subprocess.run(command, input=f">primer\n{primer_seq}", capture_output=True, text=True, check=True)
        
        hit_locations = []
        if process.stdout.strip():
            for line in process.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    hit_locations.append((parts[0], int(parts[1]), int(parts[2])))
        
        hit_count = len(hit_locations)
        is_specific = hit_count <= strict_threshold
        
        warning_message = None
        if not is_specific:
            locations_str = ', '.join([f"{loc[0]}:{loc[1]}-{loc[2]}" for loc in hit_locations[:3]])
            if hit_count > 3:
                locations_str += f" (+{hit_count - 3} more)"
            warning_message = f"Non-specific primer ({hit_count} hits): {locations_str}"
        
        return hit_count, hit_locations, is_specific, warning_message
        
    except subprocess.CalledProcessError as e:
        # BLAST failed, return safe defaults
        return 0, [], True, f"BLAST check failed: {e}"

def write_design_output_files(all_csv_rows, all_bed_rows, output_prefix, final_warnings, failed_targets_initial):
    if not all_csv_rows:
        print(f"\nNo specific primers were successfully designed for '{output_prefix}'.")
        return
        
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    csv_file, bed_file = f"{output_prefix}.csv", f"{output_prefix}.bed"
    trimming_fasta = f"{output_prefix}_trimming.fasta"
    
    csv_headers = [
        'target_id', 'pair_rank', 'flags',
        'fwd_primer_tailed', 'rev_primer_tailed',
        'fwd_primer_seq', 'rev_primer_seq',
        'fwd_primer_tm', 'rev_primer_tm', 
        'amplicon_size', 'specificity_hits'
    ]
    
    # Write CSV
    with open(csv_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=csv_headers)
        writer.writeheader()
        writer.writerows(all_csv_rows)
    
    # Write BED
    with open(bed_file, 'w') as bedf:
        bedf.writelines(all_bed_rows)
        
    # Write Trimming FASTA (Universal Tail + Target Specific Primer)
    # This ignores the hairpin clamp as it is not present in the reads.
    with open(trimming_fasta, 'w') as fastaf:
        for row in all_csv_rows:
            fwd_name = f"{row['target_id']}_F"
            rev_name = f"{row['target_id']}_R"
            fwd_seq_full = FWD_P5_TRUNCATED + row['fwd_primer_seq']
            rev_seq_full = REV_P7_TRUNCATED + row['rev_primer_seq']
            fastaf.write(f">{fwd_name}\n{fwd_seq_full}\n")
            fastaf.write(f">{rev_name}\n{rev_seq_full}\n")

    # (v8.0 Ready) Write Thermodynamic/Dimerization Report
    dimer_report_file = f"{output_prefix}.dimers.txt"
    with open(dimer_report_file, 'w') as dimf:
        dimf.write(f"--- Dimerization & Thermodynamic Report for {output_prefix} ---\n")
        dimf.write(f"Panel Size: {len(all_csv_rows)} amplicons\n\n")
        
        if final_warnings:
            dimf.write("WARNINGS / RESIDUAL CLASHES:\n")
            # Sort warnings by dG (worst first)
            try:
                sorted_warnings = sorted(list(set(final_warnings)), key=lambda x: x[1])
                for msg, dg in sorted_warnings:
                    dimf.write(f"  [{dg:.2f} kcal/mol] {msg}\n")
            except:
                for warning in final_warnings:
                    dimf.write(f"  {warning}\n")
        else:
            dimf.write("Thermodynamic Success: No significant dimerization predicted (dG thresholds respected).\n")

    print(f"\nResults for {len(all_csv_rows)} specific targets saved to:")
    print(f"   - {csv_file}")
    print(f"   - {bed_file}")
    print(f"   - {trimming_fasta}")
    print(f"   - {dimer_report_file}")
    
    log_file = f"{output_prefix}.log.txt"
    with open(log_file, 'w') as logf:
        logf.write(f"--- Design Log for {output_prefix} ---\n\n")
        
        if final_warnings:
            logf.write("--- Final Compatibility Warnings ---\n")
            try:
                sorted_warnings = sorted(list(set(final_warnings)), key=lambda x: x[1])
            except TypeError:
                sorted_warnings = final_warnings 
            for warning_msg, dg_val in sorted_warnings:
                logf.write(f"   - {warning_msg}\n")
        
        if failed_targets_initial:
            logf.write("\n--- Targets That Failed Initial Design (All Pools) ---\n")
            for reason in sorted(list(set(failed_targets_initial))):
                logf.write(f"   - {reason}\n")
                
        if not final_warnings and not failed_targets_initial:
            logf.write("Design complete. All targets were successful and compatible.\n")
            
    print(f"A detailed log of warnings has been saved to '{log_file}'")

def process_single_target(target_id, genome_records, gene_coords, blast_db, force_multiprime, oligo_format, add_hairpin_clamp, rescue=False, num_candidates=25):
    
    # We define the "ideal" size here to use for rescue calculations
    ideal_min_product_size = 150 
    
    # Define window shifts for rescue mode
    # Motif jump is ~primer size; Regional jump is 30% of target amplicon size
    motif_jump = 25
    regional_jump = int(ideal_min_product_size * 0.3)
    shifts = [0, -motif_jump, motif_jump, -regional_jump, regional_jump] if rescue else [0]
    
    for shift in shifts:
        if shift != 0:
            print(f"  -> Rescue Attempt: Shifting window by {shift}bp for {target_id}...")
            
        target_sequence, target_info, error = extract_target_sequence(genome_records, gene_coords, target_id, offset_bp=shift)
        if not target_sequence:
            continue

        all_specific_pairs_for_target = []
        seq_len = len(target_sequence)

        def find_valid_pairs(primer_results, strategy_name, current_shift):
            num_returned = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if num_returned == 0:
                return False 

            found_at_least_one = False
            for i in range(num_returned):
                fwd_seq = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
                rev_seq = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']
                
                try:
                    fwd_hits, _, fwd_specific, fwd_warning = run_blast_specificity_check(fwd_seq, blast_db)
                    rev_hits, _, rev_specific, rev_warning = run_blast_specificity_check(rev_seq, blast_db)
                except Exception as e:
                    continue 

                specificity_ok = (fwd_specific and rev_specific)
                if not force_multiprime and not specificity_ok:
                    continue 
                
                # 1. Check self-dimerization (full length)
                if get_homodimer_dg(fwd_seq) < MAX_SELF_DIMER_DG or get_homodimer_dg(rev_seq) < MAX_SELF_DIMER_DG:
                    continue

                # 2. Check internal heterodimer (full length between F and R)
                if get_heterodimer_dg(fwd_seq, rev_seq) < MAX_HETERODIMER_DG:
                    continue

                # 3. Check 3' end stability (The "Extender" check)
                dg_fwd_rev = get_end_stability_dg(fwd_seq, rev_seq)
                dg_rev_fwd = get_end_stability_dg(rev_seq, fwd_seq)
                if min(dg_fwd_rev, dg_rev_fwd) < END_STABILITY_DG_THRESHOLD:
                    continue 

                found_at_least_one = True 
                fwd_tm = primer_results[f'PRIMER_LEFT_{i}_TM']
                rev_tm = primer_results[f'PRIMER_RIGHT_{i}_TM']
                flags = []
                if current_shift != 0: flags.append(f"Rescued_{current_shift}bp")
                
                if oligo_format == 'fwd_tailed':
                    fwd_primer_tailed = FWD_P5_TRUNCATED + fwd_seq
                    rev_primer_tailed = REV_P7_TRUNCATED + rev_seq
                    if add_hairpin_clamp:
                        flags.append("HairpinClamp")
                        fwd_primer_tailed = add_iterative_hairpin_clamp(fwd_primer_tailed, fwd_seq, HAIRPIN_STEM_TARGET_DG, HAIRPIN_CLAMP_MAX_LEN)
                        rev_primer_tailed = add_iterative_hairpin_clamp(rev_primer_tailed, rev_seq, HAIRPIN_STEM_TARGET_DG, HAIRPIN_CLAMP_MAX_LEN)
                else:
                    fwd_rc = str(Seq(fwd_seq).reverse_complement())
                    rev_rc = str(Seq(rev_seq).reverse_complement())
                    fwd_primer_tailed = fwd_rc + FWD_RC_P5_TRUNCATED
                    rev_primer_tailed = rev_rc + REV_RC_P7_TRUNCATED

                pair_rank_str = f"{i} ({strategy_name})"
                csv_row = {
                    'target_id': target_id, 'pair_rank': pair_rank_str, 
                    'flags': ";".join(flags) if flags else "OK", 
                    'fwd_primer_tailed': fwd_primer_tailed, 'rev_primer_tailed': rev_primer_tailed,
                    'fwd_primer_seq': fwd_seq, 'rev_primer_seq': rev_seq,
                    'fwd_primer_tm': f"{fwd_tm:.2f}", 'rev_primer_tm': f"{rev_tm:.2f}",
                    'amplicon_size': primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                    'specificity_hits': f"F:{fwd_hits}, R:{rev_hits}"
                }
                try:
                    ideal_strategy_settings = {} 
                    # ideal_min_product_size is now defined at the top of the function
                    product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                    fwd_start_in_slice = primer_results[f'PRIMER_LEFT_{i}'][0]
                    amp_start = target_info['start'] - 1 + current_shift + fwd_start_in_slice
                    amp_end = amp_start + product_size
                    bed_row = f"{target_info['contig']}\t{amp_start}\t{amp_end}\t{target_id}_amplicon_{pair_rank_str}\t0\t+\n"
                    all_specific_pairs_for_target.append({'csv_row': csv_row, 'bed_row': bed_row})
                except Exception as e:
                    print(f"Error processing primer pair {i} for {target_id}: {e}")
                    continue
            
            return found_at_least_one

        try:
            primer_results = design_primers_for_sequence(target_sequence, target_id, {}, num_candidates=num_candidates)
            find_valid_pairs(primer_results, "Strategy 0 (Ideal)", shift)
        except: pass

        if not all_specific_pairs_for_target:
            fallback_strategies = [
                {'name': 'Strategy 1 (Longer)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[250, 350]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}},
                {'name': 'Strategy 2 (Shorter)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[100, 150]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}},
                {'name': 'Strategy 3 (Relaxed Tm)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': 55.0, 'PRIMER_MAX_TM': 65.0}},
            ]
            for strategy in fallback_strategies:
                try:
                    primer_results = design_primers_for_sequence(target_sequence, target_id, strategy['settings'], num_candidates=num_candidates)
                    if find_valid_pairs(primer_results, strategy['name'], shift): break
                except: continue

        if all_specific_pairs_for_target:
            return (target_id, all_specific_pairs_for_target, None)

    return (target_id, [], f"Could not find a specific primer pair for {target_id} even after rescue shifts.")

def is_compatible(candidate_pair, existing_pool_data, check_overlap):
    """Checks if a candidate pair is compatible with all primers and amplicons in the pool."""
    c_fwd = candidate_pair['csv_row']['fwd_primer_seq']
    c_rev = candidate_pair['csv_row']['rev_primer_seq']
    
    # 1. Internal checks for the candidate itself
    if min(get_homodimer_dg(c_fwd), get_homodimer_dg(c_rev)) < MAX_SELF_DIMER_DG:
        return False, "Self-dimer"
    if get_heterodimer_dg(c_fwd, c_rev) < MAX_HETERODIMER_DG:
        return False, "Internal FR-dimer"
    
    # 2. Genomic Overlap Check (Crucial for Tiled Designs)
    if check_overlap:
        c_bed = candidate_pair['bed_row'].split('\t')
        c_contig, c_start, c_end = c_bed[0], int(c_bed[1]), int(c_bed[2])
        for p_data in existing_pool_data:
            p_bed = p_data['bed_row'].split('\t')
            p_contig, p_start, p_end = p_bed[0], int(p_bed[1]), int(p_bed[2])
            if c_contig == p_contig:
                # Standard overlap check: (StartA < EndB) and (EndA > StartB)
                if c_start < p_end and c_end > p_start:
                    return False, "Genomic overlap"

    # 3. Cross-thermodynamic checks against the pool
    for p_data in existing_pool_data:
        p_seqs = [p_data['csv_row']['fwd_primer_seq'], p_data['csv_row']['rev_primer_seq']]
        for p_seq in p_seqs:
            for c_seq in [c_fwd, c_rev]:
                # 3' End Stability (Bi-directional)
                if get_end_stability_dg(c_seq, p_seq) < END_STABILITY_DG_THRESHOLD or \
                   get_end_stability_dg(p_seq, c_seq) < END_STABILITY_DG_THRESHOLD:
                    return False, "3' End Dimer"
                
                # Full Heterodimer
                if get_heterodimer_dg(c_seq, p_seq) < MAX_HETERODIMER_DG:
                    return False, "Internal Heterodimer"
                
    return True, None

def pack_targets_first_fit(target_to_primers_map, overlap_detected, monte_carlo_iterations=20):
    """
    (v8.5 Ready) Monte Carlo Packing Engine.
    Runs First-Fit multiple times with randomized target orderings to find the global minimum.
    """
    import random
    
    def get_target_coords(tid):
        bed = target_to_primers_map[tid][0]['bed_row'].split('\t')
        return (bed[0], int(bed[1]))

    all_tids = list(target_to_primers_map.keys())
    best_pools = None
    min_pool_count = float('inf')

    print(f"\n--- Packing {len(all_tids)} targets (Monte Carlo Iterations: {monte_carlo_iterations}) ---")
    
    for i in range(monte_carlo_iterations):
        # On iteration 0, use genomic order (best heuristic). Subsequent iterations use random.
        current_shuffled_tids = sorted(all_tids, key=get_target_coords) if i == 0 else random.sample(all_tids, len(all_tids))
        
        current_pools = [] # List of lists of primer_data
        
        for tid in current_shuffled_tids:
            candidates = target_to_primers_map[tid]
            placed = False
            
            # Optimization: Try existing pools
            for pool_set in current_pools:
                for rank_idx, cand in enumerate(candidates):
                    compatible, _ = is_compatible(cand, pool_set, overlap_detected)
                    if compatible:
                        # (v8.6 Ready) Track the rank actually used in the packing result
                        cand_with_rank = cand.copy()
                        cand_with_rank['rank_used'] = rank_idx
                        pool_set.append(cand_with_rank)
                        placed = True
                        break
                if placed: break
            
            if not placed:
                # New pool
                best_cand = candidates[0].copy()
                best_cand['rank_used'] = 0
                current_pools.append([best_cand])
        
        if len(current_pools) < min_pool_count:
            min_pool_count = len(current_pools)
            best_pools = [p[:] for p in current_pools]
            if min_pool_count <= (2 if overlap_detected else 1): 
                break # Found the theoretical minimum
            
        if i % 5 == 0 and i > 0:
            print(f"   -> Iter {i}: Min pools found so far: {min_pool_count}")

    return best_pools

def check_amplicon_overlap(best_primer_pairs):
    amplicons_by_contig = {}
    
    for pair_data in best_primer_pairs:
        bed_row = pair_data['bed_row'].strip().split('\t')
        contig, start, end = bed_row[0], int(bed_row[1]), int(bed_row[2])
        
        if contig not in amplicons_by_contig:
            amplicons_by_contig[contig] = []
        amplicons_by_contig[contig].append((start, end))
        
    for contig, coords in amplicons_by_contig.items():
        if len(coords) < 2:
            continue
        
        sorted_coords = sorted(coords, key=lambda x: x[0])
        
        for i in range(len(sorted_coords) - 1):
            current_end = sorted_coords[i][1]
            next_start = sorted_coords[i+1][0]
            
            if current_end > next_start: # Overlap!
                return True 
                
    return False 

def run_design_mode(args):
    print("Running in 'Full Design' mode...")
    failed_targets_initial = []

    try:
        create_blast_db_if_needed(args.genome, args.blast_db)
        base, ext = os.path.splitext(args.genome)
        if base.endswith('.cleaned'):
            cleaned_genome_path = args.genome
        else:
            cleaned_genome_path = f"{base}.cleaned{ext}"
        print(f"Loading CLEANED genome from '{cleaned_genome_path}'...")
        genome_records = SeqIO.to_dict(SeqIO.parse(cleaned_genome_path, "fasta"))
        print(f"Parsing GFF file '{args.gff}'...")
        gene_coords = parse_gff(args.gff)
        print(f"Reading target IDs from '{args.target_file}'...")
        target_ids = read_lines_from_file(args.target_file)

        process_func = partial(process_single_target,
                               genome_records=genome_records,
                               gene_coords=gene_coords,
                               blast_db=args.blast_db,
                               force_multiprime=args.force_multiprime,
                               oligo_format=args.oligo_format, 
                               add_hairpin_clamp=args.add_hairpin_clamp,
                               rescue=args.rescue,
                               num_candidates=args.num_candidates) 

        num_workers = os.cpu_count()
        print(f"Starting parallel processing with {num_workers} workers for {len(target_ids)} targets...")
        
        results = []
        with multiprocessing.Pool(processes=num_workers) as pool:
            results = list(tqdm(pool.imap_unordered(process_func, target_ids), total=len(target_ids), desc="Designing Primers"))

        target_to_primers_map = {}
        for target_id, primer_list, error in results:
            if primer_list:
                sorted_list = sorted(primer_list, key=lambda x: x['csv_row']['pair_rank'])
                target_to_primers_map[target_id] = sorted_list
            else:
                if error: 
                    failed_targets_initial.append(error)
        
        if not target_to_primers_map:
             print("\nNo primers were found for any target with the specified criteria. Exiting.") 
             return
             
        best_primer_pairs = []
        # (v8.0 Ready) Iterate over keys, not original list
        for target_id in target_to_primers_map:
             best_primer_pairs.append(target_to_primers_map[target_id][0]) 
        
        overlap_detected = check_amplicon_overlap(best_primer_pairs)
        
        if overlap_detected:
            print("\nOverlap detected between amplicons. This is a 'tiled' design.")
            print("Forcing interleaved multi-pool design to prevent on-target cross-priming.")
            run_minimum_pool_logic = True
        
        else:
            print("\nAmplicons are 'sparse' (non-overlapping).")
            if args.force_single_pool:
                print("User specified '--force-single-pool'. Attempting to find best-available single set.")
                run_minimum_pool_logic = False
            else:
                print("Defaulting to safest method: finding minimum number of compatible pools.")
                run_minimum_pool_logic = True
        
        if run_minimum_pool_logic:
            # (v8.5 Ready) Pass monte_carlo_iterations and force_sparse
            check_overlap = overlap_detected or args.force_sparse
            final_pools = pack_targets_first_fit(target_to_primers_map, check_overlap, args.monte_carlo_iters)
            
            print(f"\nSUCCESS: Compressed {len(target_to_primers_map)} targets into {len(final_pools)} pools.")
            
            # (v8.6 Ready) Report candidate usage statistics
            all_ranks = [p['rank_used'] for pool in final_pools for p in pool]
            avg_rank = sum(all_ranks) / len(all_ranks)
            max_rank = max(all_ranks)
            print(f"Candidate usage: Avg Rank={avg_rank:.1f}, Max Rank={max_rank}")

            for i, pool_set in enumerate(final_pools):
                pool_name = i+1 if isinstance(i, int) else i # Safety
                pool_csv = [row['csv_row'] for row in pool_set]
                pool_bed = [row['bed_row'] for row in pool_set]
                write_design_output_files(pool_csv, pool_bed, f"{args.output_prefix}_pool_{pool_name}", [], failed_targets_initial if i == 0 else [])
        
        else:
            # Forced single pool - simplified
            check_overlap = overlap_detected or args.force_sparse
            final_pools = pack_targets_first_fit(target_to_primers_map, check_overlap, args.monte_carlo_iters)
            pool_set = final_pools[0]
            write_design_output_files([r['csv_row'] for r in pool_set], [r['bed_row'] for r in pool_set], args.output_prefix, [], failed_targets_initial)

        if failed_targets_initial:
            print("\n--- Targets That Failed Initial Design (All Pools) ---")
            for reason in sorted(list(set(failed_targets_initial))):
                print(f"   - {reason}")

    except (FileNotFoundError, subprocess.CalledProcessError, RuntimeError, ValueError) as e:
        print(f"\nAn error occurred: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

# --- Mode 2: Tail-Only Pipeline Functions ---

def write_tail_only_csv(all_csv_rows, output_prefix):
    if not all_csv_rows:
        print("\nNo primers were processed.")
        return
    
    csv_file = f"{output_prefix}.csv"
    csv_headers = [
        'pair_id', 'fwd_primer_tailed', 'rev_primer_tailed',
        'fwd_primer_seq', 'rev_primer_seq'
    ]
    with open(csv_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=csv_headers)
        writer.writeheader()
        writer.writerows(all_csv_rows)
    print(f"\nSuccessfully tailed {len(all_csv_rows)} primer pairs.")
    print(f"Results saved to '{csv_file}'")

def run_tail_only_mode(args):
    print(f"Reading forward primers from: {args.tail_fwd_file}")
    print(f"Reading reverse primers from: {args.tail_rev_file}")
    
    fwd_primers = read_lines_from_file(args.tail_fwd_file)
    rev_primers = read_lines_from_file(args.tail_rev_file)
    
    if len(fwd_primers) != len(rev_primers):
        print("\nWarning: The number of forward and reverse primers does not match.")
        print(f"   Forward primers found: {len(fwd_primers)}")
        print(f"   Reverse primers found: {len(rev_primers)}")
        print("Processing the minimum number of pairs.")
    
    all_csv_rows = []
    num_pairs = min(len(fwd_primers), len(rev_primers))
    
    for i in range(num_pairs):
        fwd_seq = fwd_primers[i]
        rev_seq = rev_primers[i]
        
        if args.oligo_format == 'fwd_tailed':
            fwd_tailed = FWD_P5_TRUNCATED + fwd_seq
            rev_tailed = REV_P7_TRUNCATED + rev_seq
        else:
            fwd_rc = str(Seq(fwd_seq).reverse_complement())
            rev_rc = str(Seq(rev_seq).reverse_complement())
            fwd_tailed = fwd_rc + FWD_RC_P5_TRUNCATED
            rev_tailed = rev_rc + REV_RC_P7_TRUNCATED
        
        all_csv_rows.append({
            'pair_id': f"pair_{i+1}",
            'fwd_primer_tailed': fwd_tailed,
            'rev_primer_tailed': rev_tailed,
            'fwd_primer_seq': fwd_seq,
            'rev_primer_seq': rev_seq
        })
        
    write_tail_only_csv(all_csv_rows, args.output_prefix)

# --- Mode 3: Feedback Loop Pipeline Functions ---

def parse_feedback(feedback_file):
    """Loads and returns the feedback JSON data."""
    if not feedback_file:
        return None
    try:
        with open(feedback_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error parsing feedback file '{feedback_file}': {e}")
        return None

def generate_pool_recipe(feedback_data, output_prefix):
    """
    Generates a pool re-balancing recipe (CSV) based on recommendations.
    Assumes a base volume of 10uL for 1.0x concentration.
    """
    recipe_file = f"{output_prefix}_rebalancing_recipe.csv"
    base_vol = 10.0 # uL
    
    headers = ['Target', 'Current_Rel_Mean', 'Adjustment_Factor', 'Action', 'Recommended_Vol_uL']
    rows = []
    
    for rec in feedback_data.get('target_recommendations', []):
        adj = rec.get('recommended_adjustment', 1.0)
        rows.append({
            'Target': rec['target_id'],
            'Current_Rel_Mean': f"{rec['rel_mean']:.2f}x",
            'Adjustment_Factor': f"{adj:.2f}x",
            'Action': rec['action'],
            'Recommended_Vol_uL': f"{base_vol * adj:.2f}"
        })
    
    if rows:
        with open(recipe_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Re-balancing recipe saved to '{recipe_file}'")
    else:
        print("No recommendations found in feedback data. Skipping recipe generation.")

def run_feedback_mode(args):
    """
    Handles re-balancing and re-design based on feedback.
    """
    print(f"\n--- Running in Feedback Mode ---")
    feedback_data = parse_feedback(args.feedback)
    if not feedback_data:
        return

    # 1. Always generate the re-balancing recipe if feedback is present
    generate_pool_recipe(feedback_data, args.output_prefix)

    # 2. Check for targets flagged for redesign
    redesign_targets = [
        rec['target_id'] for rec in feedback_data.get('target_recommendations', [])
        if "REDESIGN" in rec.get('action', "").upper()
    ]

    if redesign_targets:
        print(f"Found {len(redesign_targets)} targets flagged for redesign: {', '.join(redesign_targets)}")
        
        # We need the original genome and GFF to redesign
        if not all([args.genome, args.gff]):
            print("Error: Redesign requested but --genome and --gff were not provided.")
            return

        # Prepare temporary target file for redesign
        temp_target_file = f"{args.output_prefix}_redesign_targets.txt"
        with open(temp_target_file, 'w') as f:
            for t in redesign_targets:
                f.write(f"{t}\n")
        
        print("Re-running design for flagged targets with relaxed parameters...")
        
        # Create a copy of args for the design run
        design_args = argparse.Namespace(**vars(args))
        design_args.target_file = temp_target_file
        design_args.output_prefix = f"{args.output_prefix}_v8_redesign"
        # Enable rescue mode automatically for redesign if requested or default to True for redesign
        design_args.rescue = True 
        run_design_mode(design_args)
        
        print(f"\nFeedback loop processed. Redesigned primers are in '{design_args.output_prefix}*'")
    else:
        print("No targets flagged for redesign. Feedback loop complete.")

# --- Main Execution Block (Router) ---

def main():
    parser = argparse.ArgumentParser(description="Design specific, adapter-tailed primers for multiple targets.")
    
    design_group = parser.add_argument_group('Mode 1: Full Design Pipeline')
    design_group.add_argument('--genome', default='common_refs/ecoli_genome.fna', help="Path to the reference genome (default: common_refs/ecoli_genome.fna).")
    design_group.add_argument('--gff', default='common_refs/genomic.gff', help="Path to a GFF file (default: common_refs/genomic.gff).")
    design_group.add_argument('--target-file', help="Path to a text file with one target gene ID per line.")
    design_group.add_argument('--blast-db', default='ncbi_data/ecoli_db', help="Prefix for the local BLAST database (default: ncbi_data/ecoli_db).")
    design_group.add_argument('--force-single-pool', action='store_true',
                              help="Force the script to find the 'best-available' single pool (default: finds minimum N-pools).")
    design_group.add_argument('--force-sparse', action='store_true',
                              help="Enforce zero genomic overlaps even if design is naturally sparse. Ensures a single-tube solution is 'clean'.")
    design_group.add_argument('--monte-carlo-iters', type=int, default=20,
                              help="Number of randomized packing iterations to find the global minimum (default: 20).")
    design_group.add_argument('--force-multiprime', action='store_true',
                              help="Allows primers that hit multiple identical locations in the genome. Useful for multi-copy genes like 16S rRNA. (Default: strict single-hit specificity).")
    
    # (v8.0 Ready) Phase 2 Pool Compression Group
    pool_group = parser.add_argument_group('Pool Compression Optimizations')
    pool_group.add_argument('--num-candidates', type=int, default=25, help='Number of primer pairs to design per target (default: 25)')
    pool_group.add_argument('--max-heterodimer-dg', type=float, default=-8.0, help='Max heterodimer dG threshold (default: -8.0)')
    pool_group.add_argument('--max-self-dimer-dg', type=float, default=-9.0, help='Max self-dimer dG threshold (default: -9.0)')
    pool_group.add_argument('--max-end-stability-dg', type=float, default=-9.0, help="Max 3' end stability dG (default: -9.0)")
    
    # (v8.0 Ready) argparse default now comes from the top-level constant
    design_group.add_argument('--max-compatibility-iterations', type=int, default=MAX_COMPATIBILITY_ITERATIONS,
                              help=f"Maximum iterations for the compatibility auto-healing algorithm. (Default: {MAX_COMPATIBILITY_ITERATIONS})")
    
    design_group.add_argument('--oligo-format', choices=['rc_tailed', 'fwd_tailed'], default='rc_tailed',
                              help="Format for 'fwd_primer_tailed' and 'rev_primer_tailed' columns. 'rc_tailed' (default) is reverse-complement + full tail. 'fwd_tailed' is forward-seq + truncated P5/P7 tail (synthesis-ready).")
    
    design_group.add_argument('--add-hairpin-clamp', action='store_true',
                              help=f"Iteratively adds a 5' clamp to 'fwd_tailed' oligos to form a hairpin with target dG ~{HAIRPIN_STEM_TARGET_DG} (stem only). REQUIRES --oligo-format fwd_tailed.")

    feedback_group = parser.add_argument_group('Mode 3: v8.0 Feedback Loop')
    feedback_group.add_argument('--feedback', help="Path to primer_balancing_feedback.json from analysis pipeline.")
    feedback_group.add_argument('--rescue', action='store_true', help="Automated Primer Rescue: Attempt to shift target window if original design fails/has zero coverage.")
    
    tail_group = parser.add_argument_group('Mode 2: Tail-Only Utility')
    tail_group.add_argument('--tail-fwd-file', help="Path to a text file with one forward primer per line.")
    tail_group.add_argument('--tail-rev-file', help="Path to a text file with one reverse primer per line.")
    
    parser.add_argument('--output-prefix', default='output/final_primers', help="Prefix for output files (default: output/final_primers).")
    
    args = parser.parse_args()

    # --- (v8.0 Ready) Globalization of Constants from Arguments ---
    global MAX_HETERODIMER_DG, MAX_SELF_DIMER_DG, END_STABILITY_DG_THRESHOLD
    MAX_HETERODIMER_DG = args.max_heterodimer_dg
    MAX_SELF_DIMER_DG = args.max_self_dimer_dg
    END_STABILITY_DG_THRESHOLD = args.max_end_stability_dg

    if args.add_hairpin_clamp and args.oligo_format != 'fwd_tailed':
        parser.error("--add-hairpin-clamp can only be used with --oligo-format fwd_tailed.")

    is_full_design_mode = all([args.genome, args.gff, args.target_file, args.blast_db])
    is_tail_only_mode = all([args.tail_fwd_file, args.tail_rev_file])
    is_feedback_mode = bool(args.feedback)

    if is_feedback_mode:
        run_feedback_mode(args)
    elif is_full_design_mode and not is_tail_only_mode:
        run_design_mode(args)
    elif is_tail_only_mode and not is_full_design_mode:
        run_tail_only_mode(args)
    elif is_full_design_mode and is_tail_only_mode:
        print("Error: Ambiguous command. Arguments for both 'Full Design' and 'Tail-Only' modes were provided.")
        print("Please choose only one mode by providing its respective arguments.")
        parser.print_help()
    else:
        if args.tail_fwd_file or args.tail_rev_file:
             print("Error: For 'Tail-Only' mode, you MUST provide BOTH --tail-fwd-file and --tail-rev-file.")
        elif args.feedback:
             pass # Handled above
        else:
            print("Error: You must provide the correct arguments for a mode.")
            print("\nFor 'Full Design' mode, you MUST provide:")
            print("   --genome, --gff, --target-file, and --blast-db")
            print("\nFor 'Tail-Only' mode, you MUST provide:")
            print("   --tail-fwd-file and --tail-rev-file")
            print("\nFor 'v8.0 Feedback' mode, you MUST provide:")
            print("   --feedback [json_file]")
        parser.print_help()

if __name__ == "__main__":
    main()
