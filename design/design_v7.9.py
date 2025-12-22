import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq # Added for reverse complement
import primer3
import os
import csv
import subprocess
import shutil
import multiprocessing
import itertools
from functools import partial
from tqdm import tqdm
import re # Import re for parsing

# --- (v7.7) Refactored Hardcoded Adapter Tails ---
# Tails are the full length from index to insert
#FWD_P5_TRUNCATED = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT' 33-mer
#REV_P7_TRUNCATED = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' 34-mer

# --- (v7.9) Refactored Hardcoded Adapter Tails ---
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
# (v7.9) Created a single source of truth for min primer size
IDEAL_PRIMER_MIN_SIZE = 19 

# --- Large-Panel Logic ---
MAX_CLASH_RECOMMENDATION = 5
# (v7.9) Set the default for iterations to 5000, as discussed
MAX_COMPATIBILITY_ITERATIONS = 5000

# --- (v7.9) Constants for Hairpin Clamp Logic ---
# Calibrated based on v7.5.1 results.
HAIRPIN_STEM_TARGET_DG = -11.9 # (kcal/mol) Calibrated target dG for the *stem* interaction
HAIRPIN_CLAMP_MAX_LEN = IDEAL_PRIMER_MIN_SIZE  # (v7.9) Linked to min primer size

# --- Core Helper Functions ---

def read_lines_from_file(file_path):
    """Reads a list of items from a file, one per line, skipping empty lines."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# --- (v7.8) Corrected Helper Function for Hairpin Clamps ---

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

def extract_target_sequence(genome_records, gene_coords, target_id):
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
    start, end = target_info['start'] - 1, target_info['end']
    target_seq = contig_seq[start:end]
    return str(target_seq), target_info, None

def design_primers_for_sequence(sequence, target_id, strategy_settings):
    seq_args = {'SEQUENCE_ID': target_id, 'SEQUENCE_TEMPLATE': sequence}
    
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': IDEAL_PRIMER_MIN_SIZE, # (v7.9) Use constant
        'PRIMER_MAX_SIZE': 21,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': IDEAL_TM_MIN,
        'PRIMER_MAX_TM': IDEAL_TM_MAX,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]],
        'PRIMER_NUM_RETURN': 20
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
    with open(trimming_fasta, 'w') as tf:
        for row in all_csv_rows:
            target_id = row['target_id']
            fwd_seq = row['fwd_primer_seq']
            rev_seq = row['rev_primer_seq']
            
            # Universal Tail + Target Specific Primer
            tf.write(f">{target_id}_F\n{FWD_P5_TRUNCATED}{fwd_seq}\n")
            tf.write(f">{target_id}_R\n{REV_P7_TRUNCATED}{rev_seq}\n")

    print(f"\nResults for {len(all_csv_rows)} specific targets saved to:")
    print(f"   - {csv_file}")
    print(f"   - {bed_file}")
    print(f"   - {trimming_fasta}")
    
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

def process_single_target(target_id, genome_records, gene_coords, blast_db, force_multiprime, oligo_format, add_hairpin_clamp):
    
    target_sequence, target_info, error = extract_target_sequence(genome_records, gene_coords, target_id)
    if not target_sequence:
        return (target_id, [], error) 

    all_specific_pairs_for_target = []
    seq_len = len(target_sequence)

    def find_valid_pairs(primer_results, strategy_name):
        num_returned = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if num_returned == 0:
            return False 

        found_at_least_one = False
        for i in range(num_returned):
            fwd_seq = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_seq = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']
            
            try:
                fwd_hits = run_blast_specificity_check(fwd_seq, blast_db)
                rev_hits = run_blast_specificity_check(rev_seq, blast_db)
            except Exception as e:
                print(f"Warning: BLAST check failed for {target_id} pair {i}: {e}. Skipping pair.")
                continue 

            specificity_ok = (fwd_hits == 1 and rev_hits == 1)
            
            if not force_multiprime and not specificity_ok:
                continue 
            
            fwd_self_dg = primer3.calc_end_stability(fwd_seq, fwd_seq).dg / 1000.0
            rev_self_dg = primer3.calc_end_stability(rev_seq, rev_seq).dg / 1000.0
            fwd_rev_dg1 = primer3.calc_end_stability(fwd_seq, rev_seq).dg / 1000.0
            fwd_rev_dg2 = primer3.calc_end_stability(rev_seq, fwd_seq).dg / 1000.0
            
            worst_self_interaction_dg = min(fwd_self_dg, rev_self_dg, fwd_rev_dg1, fwd_rev_dg2)

            if worst_self_interaction_dg < END_STABILITY_DG_THRESHOLD:
                continue 

            found_at_least_one = True 
            fwd_tm = primer_results[f'PRIMER_LEFT_{i}_TM']
            rev_tm = primer_results[f'PRIMER_RIGHT_{i}_TM']
            flags = []
            
            if oligo_format == 'fwd_tailed':
                fwd_primer_tailed = FWD_P5_TRUNCATED + fwd_seq
                rev_primer_tailed = REV_P7_TRUNCATED + rev_seq
                
                if add_hairpin_clamp:
                    flags.append("HairpinClamp")
                    fwd_primer_tailed = add_iterative_hairpin_clamp(
                        fwd_primer_tailed, fwd_seq, HAIRPIN_STEM_TARGET_DG, HAIRPIN_CLAMP_MAX_LEN
                    )
                    rev_primer_tailed = add_iterative_hairpin_clamp(
                        rev_primer_tailed, rev_seq, HAIRPIN_STEM_TARGET_DG, HAIRPIN_CLAMP_MAX_LEN
                    )
            else:
                fwd_rc = str(Seq(fwd_seq).reverse_complement())
                rev_rc = str(Seq(rev_seq).reverse_complement())
                fwd_primer_tailed = fwd_rc + FWD_RC_P5_TRUNCATED
                rev_primer_tailed = rev_rc + REV_RC_P7_TRUNCATED

            if fwd_tm < IDEAL_TM_MIN or rev_tm < IDEAL_TM_MIN:
                flags.append("Low_Tm")
            if fwd_tm > IDEAL_TM_MAX or rev_tm > IDEAL_TM_MAX:
                flags.append("High_Tm")
            
            if force_multiprime and not specificity_ok:
                flags.append(f"Multi_Hit_FWD:{fwd_hits}_REV:{rev_hits}")

            pair_rank_str = f"{i} ({strategy_name})"

            csv_row = {
                'target_id': target_id, 'pair_rank': pair_rank_str, 
                'flags': ";".join(flags) if flags else "OK", 
                'fwd_primer_tailed': fwd_primer_tailed,
                'rev_primer_tailed': rev_primer_tailed,
                'fwd_primer_seq': fwd_seq, 'rev_primer_seq': rev_seq,
                'fwd_primer_tm': f"{fwd_tm:.2f}",
                'rev_primer_tm': f"{rev_tm:.2f}",
                'amplicon_size': primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                'specificity_hits': f"F:{fwd_hits}, R:{rev_hits}"
            }
            
            product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
            fwd_start_in_gene = primer_results[f'PRIMER_LEFT_{i}'][0]
            amp_start = target_info['start'] - 1 + fwd_start_in_gene
            amp_end = amp_start + product_size
            bed_row = f"{target_info['contig']}\t{amp_start}\t{amp_end}\t{target_id}_amplicon_{pair_rank_str}\t0\t+\n"
            
            all_specific_pairs_for_target.append({'csv_row': csv_row, 'bed_row': bed_row})
        
        return found_at_least_one

    try:
        ideal_strategy_settings = {} 
        ideal_min_product_size = 150 
        
        ideal_strategy_possible = (seq_len >= ideal_min_product_size)
        
        if ideal_strategy_possible:
            primer_results = design_primers_for_sequence(target_sequence, target_id, ideal_strategy_settings)
            find_valid_pairs(primer_results, "Strategy 0 (Ideal)")

    except Exception as e:
        return (target_id, [], f"Error processing {target_id} (Strategy 0): {e}")

    if not all_specific_pairs_for_target:
        fallback_strategies = [
            {'name': 'Strategy 1 (Longer)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[250, 350]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}},
            {'name': 'Strategy 2 (Shorter)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[100, 150]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}},
            {'name': 'Strategy 3 (Relaxed Tm)', 'settings': {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': 55.0, 'PRIMER_MAX_TM': 65.0}},
        ]

        for strategy in fallback_strategies:
            strategy_name = strategy['name']
            strategy_settings = strategy['settings']
            
            product_range_list = strategy_settings.get('PRIMER_PRODUCT_SIZE_RANGE', [[150,250]])
            min_product_size = product_range_list[0][0]
            if seq_len < min_product_size:
                continue 

            try:
                primer_results = design_primers_for_sequence(target_sequence, target_id, strategy_settings)
                if find_valid_pairs(primer_results, strategy_name):
                    break 
            except Exception as e:
                print(f"Warning: Strategy {strategy_name} for {target_id} failed: {e}")
                continue 

    if not all_specific_pairs_for_target:
        return (target_id, [], f"Could not find a specific, non-hairpin primer pair for {target_id} after all strategies.")
    
    return (target_id, all_specific_pairs_for_target, None) 

def find_compatible_set(target_to_primers_map, pool_name, max_iterations):
    print(f"\n--- Starting Compatibility Check & Auto-Healing for {pool_name} ---")
    
    if len(target_to_primers_map) < 2:
        print(f"Not enough targets in {pool_name} for a multiplex check. Using best primers.")
        final_set = [primers[0] for primers in target_to_primers_map.values() if primers]
        return final_set, []

    current_indices = {target: 0 for target in target_to_primers_map}
    best_set_indices = current_indices.copy()
    min_clashes_found = float('inf')
    best_clash_list = []
    
    for iteration in range(max_iterations):
        current_primer_set = [] 
        valid_set = True
        
        for target_id, index in current_indices.items():
            if index >= len(target_to_primers_map[target_id]):
                valid_set = False
                break
            
            primer_data = target_to_primers_map[target_id][index]
            csv_row = primer_data['csv_row']
            fwd_name = f"{csv_row['target_id']}_F_{csv_row['pair_rank']}"
            rev_name = f"{csv_row['target_id']}_R_{csv_row['pair_rank']}"
            
            current_primer_set.append((fwd_name, csv_row['fwd_primer_seq'], target_id))
            current_primer_set.append((rev_name, csv_row['rev_primer_seq'], target_id))
        
        if not valid_set:
            print(f"Stopping search for {pool_name} (ran out of alternatives).")
            final_set = [target_to_primers_map[t][i] for t, i in best_set_indices.items() if i < len(target_to_primers_map[t])]
            return final_set, best_clash_list

        clashes = {} 
        all_clashes = [] 

        for (p1_name, p1_seq, p1_target), (p2_name, p2_seq, p2_target) in itertools.combinations(current_primer_set, 2):
            if p1_target == p2_target:
                continue
            try:
                dg1_on_2 = primer3.calc_end_stability(p1_seq, p2_seq).dg / 1000.0
                dg2_on_1 = primer3.calc_end_stability(p2_seq, p1_seq).dg / 1000.0
                dg_val = min(dg1_on_2, dg2_on_1)
                
                if dg_val < END_STABILITY_DG_THRESHOLD:
                    error_msg = (f"Potential 3' cross-dimer between {p1_name} and {p2_name} (dG: {dg_val:.2f})")
                    all_clashes.append((error_msg, dg_val))
                    clashes[p1_target] = clashes.get(p1_target, 0) + 1
                    clashes[p2_target] = clashes.get(p2_target, 0) + 1
            except Exception as e:
                print(f"Warning: Could not calculate 3' end stability for {p1_name}/{p2_name}. Error: {e}")

        num_clashes = len(all_clashes)
        
        if num_clashes < min_clashes_found:
            min_clashes_found = num_clashes
            best_set_indices = current_indices.copy()
            best_clash_list = all_clashes 

        if not clashes:
            print(f"SUCCESS: Found a compatible set for {pool_name} in {iteration + 1} iterations.")
            final_set = [target_to_primers_map[t][i] for t, i in current_indices.items()]
            return final_set, []
        
        worst_target = max(clashes, key=clashes.get)
        if pool_name == "Single_Pool": 
            print(f"   -> {pool_name} Iteration {iteration+1}: Found {num_clashes} clashes. Worst offender: {worst_target} ({clashes[worst_target]} clashes).")
        
        current_indices[worst_target] += 1
        
    print(f"Failed to find a perfect set for {pool_name} after {max_iterations} iterations.")
    
    worst_dg = 0.0
    if best_clash_list: 
        try:
            worst_dg = min(dg for msg, dg in best_clash_list)
        except (ValueError, TypeError):
            worst_dg = 0.0 
            
    if pool_name == "Single_Pool" or min_clashes_found > 0:
        print(f"Returning best available set for {pool_name} with {min_clashes_found} potential clashes (worst 3' dG: {worst_dg:.2f} kcal/mol).")
        
    final_set = [target_to_primers_map[t][i] for t, i in best_set_indices.items() if i < len(target_to_primers_map[t])]
    return final_set, best_clash_list 

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
                               add_hairpin_clamp=args.add_hairpin_clamp) 

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
        # (v7.9) Iterate over keys, not original list
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
            print("\nDesign Strategy: Searching for the minimum number of compatible pools...")
            
            # (v7.9) Iterate over keys, not original list
            successful_targets = sorted(list(target_to_primers_map.keys()))
            for num_pools in range(1, len(successful_targets) + 1):
                print(f"\n--- Testing N = {num_pools} {'pool' if num_pools == 1 else 'pools'} ---")
                
                pools = [{} for _ in range(num_pools)]
                for i, target_id in enumerate(successful_targets):
                    pool_index = i % num_pools
                    pools[pool_index][target_id] = target_to_primers_map[target_id]
                
                print(f"Splitting into {num_pools} pools with {[len(p) for p in pools]} targets each.")

                all_pools_succeeded = True
                all_pools_data = [] 
                
                for i, pool_targets in enumerate(pools):
                    pool_name = f"pool_{i+1}"
                    if not pool_targets:
                        continue
                    
                    pool_set, pool_clash_data = find_compatible_set(pool_targets, pool_name, max_iterations=args.max_compatibility_iterations) 
                    
                    if pool_clash_data: 
                        all_pools_succeeded = False
                        worst_dg = 0.0
                        try:
                            worst_dg = min(dg for msg, dg in pool_clash_data)
                        except (ValueError, TypeError): pass
                        print(f"Failed to find a perfect 0-clash set for {pool_name} (found {len(pool_clash_data)} clashes, worst 3' dG: {worst_dg:.2f} kcal/mol).")
                        print("Trying with N+1 pools...")
                        break 
                    
                    all_pools_data.append((pool_set, pool_clash_data, pool_name))
                
                if all_pools_succeeded:
                    print(f"\nSUCCESS: Found a perfect solution with {num_pools} pools.")
                    for pool_set, pool_clash_data, pool_name in all_pools_data:
                        pool_csv = [row['csv_row'] for row in pool_set]
                        pool_bed = [row['bed_row'] for row in pool_set]
                        write_design_output_files(pool_csv, pool_bed, f"{args.output_prefix}_{pool_name}", pool_clash_data, failed_targets_initial if pool_name == "pool_1" else [])
                    break 
            
            if not all_pools_succeeded:
                  print("\nFATAL: Could not find a perfect N-pool solution, even with N=len(targets).")

        else:
            print("\nDesign Strategy: Attempting to find a single compatible pool (forced)...")
            final_compatible_set, final_clash_data = find_compatible_set(target_to_primers_map, "Single_Pool", max_iterations=args.max_compatibility_iterations)
            
            if len(final_clash_data) > MAX_CLASH_RECOMMENDATION:
                print("\n--- RECOMMENDATION ---")
                print(f"Warning: Best set found has {len(final_clash_data)} clashes (threshold is {MAX_CLASH_RECOMMENDATION}).")
                print("This panel is at high risk of failure.")
                
            if final_compatible_set:
                final_csv_rows = [row['csv_row'] for row in final_compatible_set]
                final_bed_rows = [row['bed_row'] for row in final_compatible_set]
                write_design_output_files(final_csv_rows, final_bed_rows, args.output_prefix, final_clash_data, failed_targets_initial)
            else:
                print("\nCould not determine a final compatible set. No files will be written.")

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

# --- Main Execution Block (Router) ---

def main():
    parser = argparse.ArgumentParser(description="Design specific, adapter-tailed primers for multiple targets.")
    
    design_group = parser.add_argument_group('Mode 1: Full Design Pipeline')
    design_group.add_argument('--genome', default='../common_refs/ecoli_genome.fna', help="Path to the reference genome (default: ../common_refs/ecoli_genome.fna).")
    design_group.add_argument('--gff', default='../common_refs/genomic.gff', help="Path to a GFF file (default: ../common_refs/genomic.gff).")
    design_group.add_argument('--target-file', help="Path to a text file with one target gene ID per line.")
    design_group.add_argument('--blast-db', default='../ncbi_data/ecoli_db', help="Prefix for the local BLAST database (default: ../ncbi_data/ecoli_db).")
    design_group.add_argument('--force-single-pool', action='store_true',
                              help="Force the script to find the 'best-available' single pool (default: finds minimum N-pools).")
    design_group.add_argument('--force-multiprime', action='store_true',
                              help="Allows primers that hit multiple identical locations in the genome. Useful for multi-copy genes like 16S rRNA. (Default: strict single-hit specificity).")
    
    # (v7.9) argparse default now comes from the top-level constant
    design_group.add_argument('--max-compatibility-iterations', type=int, default=MAX_COMPATIBILITY_ITERATIONS,
                              help=f"Maximum iterations for the compatibility auto-healing algorithm. (Default: {MAX_COMPATIBILITY_ITERATIONS})")
    
    design_group.add_argument('--oligo-format', choices=['rc_tailed', 'fwd_tailed'], default='rc_tailed',
                              help="Format for 'fwd_primer_tailed' and 'rev_primer_tailed' columns. 'rc_tailed' (default) is reverse-complement + full tail. 'fwd_tailed' is forward-seq + truncated P5/P7 tail (synthesis-ready).")
    
    design_group.add_argument('--add-hairpin-clamp', action='store_true',
                              help=f"Iteratively adds a 5' clamp to 'fwd_tailed' oligos to form a hairpin with target dG ~{HAIRPIN_STEM_TARGET_DG} (stem only). REQUIRES --oligo-format fwd_tailed.")

    tail_group = parser.add_argument_group('Mode 2: Tail-Only Utility')
    tail_group.add_argument('--tail-fwd-file', help="Path to a text file with one forward primer per line.")
    tail_group.add_argument('--tail-rev-file', help="Path to a text file with one reverse primer per line.")
    
    parser.add_argument('--output-prefix', default='output/final_primers', help="Prefix for output files (default: output/final_primers).")
    
    args = parser.parse_args()

    if args.add_hairpin_clamp and args.oligo_format != 'fwd_tailed':
        parser.error("--add-hairpin-clamp can only be used with --oligo-format fwd_tailed.")

    is_full_design_mode = all([args.genome, args.gff, args.target_file, args.blast_db])
    is_tail_only_mode = all([args.tail_fwd_file, args.tail_rev_file])

    if is_full_design_mode and not is_tail_only_mode:
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
        else:
            print("Error: You must provide the correct arguments for a mode.")
            print("\nFor 'Full Design' mode, you MUST provide:")
            print("   --genome, --gff, --target-file, and --blast-db")
            print("\nFor 'Tail-Only' mode, you MUST provide:")
            print("   --tail-fwd-file and --tail-rev-file")
        parser.print_help()

if __name__ == "__main__":
    main()
