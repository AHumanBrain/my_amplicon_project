"""Multiplex PCR primer panel design and pool packing.

    python design.py --target-file genes.txt        # design + pack into pools
    python design.py --feedback coverage.json       # rebalance recipe + redesign flagged targets
"""
import argparse
import csv
import json
import logging
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
from functools import lru_cache, partial

import primer3
from Bio import SeqIO
from tqdm import tqdm

from core.biophysics import (BioConfig, add_hairpin_clamp, get_end_stability_dg,
                             get_heterodimer_dg, get_homodimer_dg)

log = logging.getLogger("design")

# Dimer thresholds (kcal/mol); overridden from the CLI and re-applied in workers.
CFG = {
    'self_dimer': BioConfig.MAX_SELF_DIMER_DG,
    'hetero_dimer': BioConfig.MAX_HETERODIMER_DG,
    'end_3p': BioConfig.END_STABILITY_DG_THRESHOLD,
}

# (index, settings). Strategy 0 is tried first; the rest are fallbacks.
STRATEGIES = [
    (0, {}),
    (1, {'PRIMER_PRODUCT_SIZE_RANGE': [[250, 350]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}),
    (2, {'PRIMER_PRODUCT_SIZE_RANGE': [[100, 150]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}),
    (3, {'PRIMER_MIN_TM': 55.0, 'PRIMER_MAX_TM': 65.0}),
]
RESCUE_SHIFTS = [0, -25, 25, -45, 45]  # bp window shifts tried with --rescue
HAIRPIN_CLAMP_MAX_LEN = 20


def configure_logging(log_file=None):
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=logging.INFO, handlers=handlers,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def sanitize_name(text):
    """BED/Picard-safe name: no whitespace or parentheses."""
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(text)).strip('_')


def parse_target_list(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip() and not ln.startswith('#')]


def target_gene_id(amplicon_name):
    """'gyrB_amplicon_S0_3' -> 'gyrB'."""
    return amplicon_name.split('_amplicon_')[0]


# --------------------------------------------------------------------------
# Inputs: genome, annotation, BLAST database
# --------------------------------------------------------------------------

def ensure_blast_db(genome_fasta, blast_db_prefix):
    if os.path.exists(f"{blast_db_prefix}.nin") or os.path.exists(f"{blast_db_prefix}.nhr"):
        return
    if not shutil.which("makeblastdb"):
        raise RuntimeError("`makeblastdb` not found; install NCBI BLAST+ (see environment.yml).")
    os.makedirs(os.path.dirname(blast_db_prefix) or '.', exist_ok=True)
    log.info("Building BLAST database '%s'", blast_db_prefix)
    try:
        subprocess.run(['makeblastdb', '-in', genome_fasta, '-dbtype', 'nucl', '-out', blast_db_prefix],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"makeblastdb failed: {e.stderr}") from e


def parse_gff(gff_file):
    """Map gene name / gene / locus_tag -> {contig, start, end, strand} (1-based, inclusive)."""
    attr_re = re.compile(r"(?:Name|gene|locus_tag)=([^;]+)")
    genes = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            for name in attr_re.findall(p[8]):
                genes.setdefault(name, {'contig': p[0], 'start': int(p[3]),
                                        'end': int(p[4]), 'strand': p[6]})
    if not genes:
        log.warning("GFF has no gene entries with Name=/gene=/locus_tag= attributes.")
    return genes


def extract_target_sequence(genome_records, gene_coords, target_id, offset_bp=0):
    """Return (sequence, gene_info, error). The window can be shifted by offset_bp."""
    if target_id not in gene_coords:
        return None, None, f"Target '{target_id}' not found in GFF."
    info = gene_coords[target_id]
    contig = info['contig']
    if contig not in genome_records:
        short = contig.split('|')[-1]
        if short not in genome_records:
            return None, None, f"Contig '{contig}' for '{target_id}' not found in FASTA."
        info = dict(info, contig=short)
    seq = genome_records[info['contig']].seq
    start = max(0, info['start'] - 1 + offset_bp)
    end = min(len(seq), info['end'] + offset_bp)
    return str(seq[start:end]), info, None


# --------------------------------------------------------------------------
# Per-target design: primer3 -> BLAST -> dimer filters
# --------------------------------------------------------------------------

def design_primers_for_sequence(sequence, target_id, strategy_settings, num_candidates=25):
    seq_args = {'SEQUENCE_ID': target_id, 'SEQUENCE_TEMPLATE': sequence}
    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 19, 'PRIMER_MAX_SIZE': 21,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': BioConfig.IDEAL_TM_MIN, 'PRIMER_MAX_TM': BioConfig.IDEAL_TM_MAX,
        'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]],
        'PRIMER_NUM_RETURN': num_candidates,
        # KAPA HiFi buffer: 2.5 mM Mg++, 1.2 mM total dNTP, SantaLucia 1998 salt correction
        'PRIMER_SALT_DIVALENT': BioConfig.SALT_DIVALENT,
        'PRIMER_DNTP_CONC': BioConfig.DNTP_CONC,
        'PRIMER_SALT_CORRECTIONS': 1,
    }
    global_args.update(strategy_settings)
    return primer3.design_primers(seq_args, global_args)


def run_blast_batch(primer_seqs, blast_db_path, strict_threshold=1, threads=1):
    """Specificity-check many primers with ONE blastn call.

    Returns {seq: (hit_count, hit_locations, is_specific, warning_or_None)}.
    Only exact, full-length matches count as hits. Fails closed on BLAST errors.
    """
    seqs = sorted(set(primer_seqs))
    if not seqs:
        return {}
    query = "".join(f">q{i}\n{q}\n" for i, q in enumerate(seqs))
    cmd = ['blastn', '-query', '-', '-db', blast_db_path, '-task', 'blastn-short',
           '-outfmt', '6 qseqid sseqid sstart send', '-perc_identity', '100',
           '-qcov_hsp_perc', '100', '-num_threads', str(threads)]
    try:
        proc = subprocess.run(cmd, input=query, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        msg = f"BLAST check failed: {e}"
        log.error(msg)
        return {q: (0, [], False, msg) for q in seqs}

    hits = {i: [] for i in range(len(seqs))}
    for line in proc.stdout.strip().splitlines():
        p = line.split('\t')
        if len(p) >= 4:
            hits[int(p[0][1:])].append((p[1], int(p[2]), int(p[3])))

    results = {}
    for i, q in enumerate(seqs):
        locs = hits[i]
        specific = len(locs) <= strict_threshold
        warning = None
        if not specific:
            shown = ', '.join(f"{c}:{s}-{e}" for c, s, e in locs[:3])
            warning = f"Non-specific primer ({len(locs)} hits): {shown}" + \
                      (f" (+{len(locs) - 3} more)" if len(locs) > 3 else "")
        results[q] = (len(locs), locs, specific, warning)
    return results


def _init_worker(cfg):
    """Pool initializer: re-apply CLI thresholds in every worker process."""
    CFG.update(cfg)


def passes_dimer_filters(fwd, rev):
    """Self-dimer, F/R heterodimer and 3' extension checks for one pair."""
    if min(get_homodimer_dg(fwd), get_homodimer_dg(rev)) < CFG['self_dimer']:
        return False
    if get_heterodimer_dg(fwd, rev) < CFG['hetero_dimer']:
        return False
    return min(get_end_stability_dg(fwd, rev), get_end_stability_dg(rev, fwd)) >= CFG['end_3p']


def process_single_target(target_id, genome_records, gene_coords, blast_db, force_multiprime,
                          hairpin_clamp=False, rescue=False, num_candidates=25):
    """Design candidate pairs for one target. Returns (target_id, candidates, error)."""
    for shift in (RESCUE_SHIFTS if rescue else [0]):
        if shift:
            log.info("Rescue: shifting window %+d bp for %s", shift, target_id)
        sequence, info, _ = extract_target_sequence(genome_records, gene_coords, target_id, shift)
        if not sequence:
            continue

        candidates, blast_cache = [], {}

        def collect(results, strategy_idx):
            n = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if n == 0:
                return
            wanted = {results[f'PRIMER_{side}_{k}_SEQUENCE'] for k in range(n) for side in ('LEFT', 'RIGHT')}
            blast_cache.update(run_blast_batch(wanted - blast_cache.keys(), blast_db))
            for i in range(n):
                fwd, rev = results[f'PRIMER_LEFT_{i}_SEQUENCE'], results[f'PRIMER_RIGHT_{i}_SEQUENCE']
                f_hits, _, f_ok, _ = blast_cache[fwd]
                r_hits, _, r_ok, _ = blast_cache[rev]
                if not force_multiprime and not (f_ok and r_ok):
                    continue
                if not passes_dimer_filters(fwd, rev):
                    continue
                flags = [f"Rescued{shift:+d}bp"] if shift else []
                fwd_oligo, rev_oligo = fwd, rev
                if hairpin_clamp:
                    flags.append("HairpinClamp")
                    fwd_oligo = add_hairpin_clamp(fwd, fwd, max_clamp_len=HAIRPIN_CLAMP_MAX_LEN)
                    rev_oligo = add_hairpin_clamp(rev, rev, max_clamp_len=HAIRPIN_CLAMP_MAX_LEN)
                start = info['start'] - 1 + shift + results[f'PRIMER_LEFT_{i}'][0]
                size = results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                candidates.append({
                    'target_id': target_id, 'strategy': strategy_idx, 'rank': i,
                    'flags': ";".join(flags) or "OK",
                    'fwd': fwd, 'rev': rev, 'fwd_oligo': fwd_oligo, 'rev_oligo': rev_oligo,
                    'fwd_tm': results[f'PRIMER_LEFT_{i}_TM'], 'rev_tm': results[f'PRIMER_RIGHT_{i}_TM'],
                    'size': size, 'contig': info['contig'], 'start': start, 'end': start + size,
                    'hits': f"F:{f_hits}, R:{r_hits}",
                    'name': f"{sanitize_name(target_id)}_amplicon_S{strategy_idx}_{i}",
                    'sort_key': (strategy_idx, i),
                })

        for idx, settings in STRATEGIES:
            try:
                collect(design_primers_for_sequence(sequence, target_id, settings, num_candidates), idx)
            except (OSError, RuntimeError) as e:  # primer3 rejects unusable templates
                log.debug("primer3 strategy %d failed for %s: %s", idx, target_id, e)
            if candidates:
                return target_id, sorted(candidates, key=lambda c: c['sort_key']), None

    return target_id, [], f"No specific primer pair found for {target_id}."


# --------------------------------------------------------------------------
# Pool packing
# --------------------------------------------------------------------------

@lru_cache(maxsize=None)
def primers_compatible(seq1, seq2, threshold_3p, threshold_full):
    """Cached cross-check of two primers: 3' extension (both ways) and full heterodimer."""
    if get_end_stability_dg(seq1, seq2) < threshold_3p or get_end_stability_dg(seq2, seq1) < threshold_3p:
        return False, "3' end dimer"
    if get_heterodimer_dg(seq1, seq2) < threshold_full:
        return False, "heterodimer"
    return True, None


def is_compatible(cand, pool, check_overlap):
    """Can `cand` join `pool`? Checks amplicon overlap (optional) and all primer cross-dimers."""
    if check_overlap:
        for p in pool:
            if p['contig'] == cand['contig'] and cand['start'] < p['end'] and cand['end'] > p['start']:
                return False, "amplicon overlap"
    for p in pool:
        for ps in (p['fwd'], p['rev']):
            for cs in (cand['fwd'], cand['rev']):
                ok, why = primers_compatible(cs, ps, CFG['end_3p'], CFG['hetero_dimer'])
                if not ok:
                    return False, why
    return True, None


def has_amplicon_overlap(candidates):
    by_contig = {}
    for c in candidates:
        by_contig.setdefault(c['contig'], []).append((c['start'], c['end']))
    for spans in by_contig.values():
        spans.sort()
        if any(spans[i][1] > spans[i + 1][0] for i in range(len(spans) - 1)):
            return True
    return False


def pack_targets_first_fit(target_to_cands, check_overlap, iterations=20, seed=0):
    """Monte Carlo first-fit: repeat with shuffled target order, keep the fewest-pool result.

    Iteration 0 uses genomic order. Each placed candidate gets 'rank_used' = its index in
    that target's ranked candidate list.
    """
    rng = random.Random(seed)
    tids = list(target_to_cands)
    genomic = sorted(tids, key=lambda t: (target_to_cands[t][0]['contig'], target_to_cands[t][0]['start']))
    best, best_n = None, float('inf')
    log.info("Packing %d targets (%d iterations)", len(tids), iterations)

    for it in range(iterations):
        order = genomic if it == 0 else rng.sample(tids, len(tids))
        pools = []
        for tid in order:
            cands = target_to_cands[tid]
            placed = False
            for pool in pools:
                for rank, cand in enumerate(cands):
                    if is_compatible(cand, pool, check_overlap)[0]:
                        pool.append(dict(cand, rank_used=rank))
                        placed = True
                        break
                if placed:
                    break
            if not placed:  # new pool: seed with the candidate that has the weakest self-dimer
                best_c = max(cands, key=lambda c: min(get_homodimer_dg(c['fwd']), get_homodimer_dg(c['rev'])))
                pools.append([dict(best_c, rank_used=cands.index(best_c))])
        if len(pools) < best_n:
            best, best_n = [p[:] for p in pools], len(pools)
            if best_n <= (2 if check_overlap else 1):
                break  # cannot do better
    return best


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

CSV_FIELDS = ['pool', 'target_id', 'name', 'flags', 'fwd_oligo', 'rev_oligo', 'fwd_primer', 'rev_primer',
              'fwd_tm', 'rev_tm', 'amplicon_size', 'contig', 'start', 'end', 'specificity_hits']


def write_outputs(pools, prefix, failed):
    """Per pool: .csv, .bed, _trimming.fasta. Overall: _order.csv (pool,name,sequence) and failures."""
    os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
    with open(f"{prefix}_order.csv", 'w', newline='') as of:
        order = csv.writer(of)
        order.writerow(['pool', 'name', 'sequence'])
        for n, pool in enumerate(pools, 1):
            rows = [{
                'pool': n, 'target_id': c['target_id'], 'name': c['name'], 'flags': c['flags'],
                'fwd_oligo': c['fwd_oligo'], 'rev_oligo': c['rev_oligo'],
                'fwd_primer': c['fwd'], 'rev_primer': c['rev'],
                'fwd_tm': f"{c['fwd_tm']:.2f}", 'rev_tm': f"{c['rev_tm']:.2f}",
                'amplicon_size': c['size'], 'contig': c['contig'], 'start': c['start'], 'end': c['end'],
                'specificity_hits': c['hits']} for c in pool]
            with open(f"{prefix}_pool_{n}.csv", 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
                w.writeheader()
                w.writerows(rows)
            with open(f"{prefix}_pool_{n}.bed", 'w') as f:
                for c in pool:
                    f.write(f"{c['contig']}\t{c['start']}\t{c['end']}\t{c['name']}\t0\t+\n")
            # Reads begin at the primer; <target>_F / <target>_R names let the pipeline split orientations.
            with open(f"{prefix}_pool_{n}_trimming.fasta", 'w') as f:
                for c in pool:
                    f.write(f">{c['target_id']}_F\n{c['fwd_oligo']}\n>{c['target_id']}_R\n{c['rev_oligo']}\n")
            for c in pool:
                order.writerow([n, f"{c['name']}_F", c['fwd_oligo']])
                order.writerow([n, f"{c['name']}_R", c['rev_oligo']])
    if failed:
        with open(f"{prefix}_failed_targets.txt", 'w') as f:
            f.write("\n".join(sorted(set(failed))) + "\n")
    log.info("Wrote %d pool(s) to %s_pool_N.* and %s_order.csv", len(pools), prefix, prefix)


# --------------------------------------------------------------------------
# Mode 1: design
# --------------------------------------------------------------------------

def run_design_mode(args):
    ensure_blast_db(args.genome, args.blast_db)
    genome_records = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    gene_coords = parse_gff(args.gff)
    target_ids = parse_target_list(args.target_file)

    worker = partial(process_single_target, genome_records=genome_records, gene_coords=gene_coords,
                     blast_db=args.blast_db, force_multiprime=args.force_multiprime,
                     hairpin_clamp=args.add_hairpin_clamp, rescue=args.rescue,
                     num_candidates=args.num_candidates)
    with multiprocessing.Pool(os.cpu_count(), initializer=_init_worker, initargs=(dict(CFG),)) as pool:
        results = list(tqdm(pool.imap_unordered(worker, target_ids), total=len(target_ids),
                            desc="Designing primers"))

    target_to_cands, failed = {}, []
    for tid, cands, err in results:
        if cands:
            target_to_cands[tid] = cands
        elif err:
            failed.append(err)
    if not target_to_cands:
        log.error("No primers found for any target.")
        sys.exit(1)

    overlap = has_amplicon_overlap([c[0] for c in target_to_cands.values()])
    log.info("Design is %s", "tiled (overlapping amplicons): packing into non-overlapping pools"
             if overlap else "sparse (non-overlapping amplicons)")
    pools = pack_targets_first_fit(target_to_cands, overlap or args.force_sparse,
                                   args.monte_carlo_iters, args.seed)

    if args.force_single_pool and len(pools) > 1:
        kept = max(pools, key=len)
        left = sorted(set(target_to_cands) - {c['target_id'] for c in kept})
        msg = f"--force-single-pool: {len(left)} target(s) could not share the tube: {', '.join(left)}"
        log.warning(msg)
        failed.append(msg)
        pools = [kept]

    ranks = [c['rank_used'] for p in pools for c in p]
    log.info("%d targets -> %d pool(s) %s; candidate rank used: mean %.1f, max %d",
             sum(len(p) for p in pools), len(pools), [len(p) for p in pools], sum(ranks) / len(ranks), max(ranks))
    write_outputs(pools, args.output_prefix, failed)
    for reason in sorted(set(failed)):
        log.warning("Failed: %s", reason)


# --------------------------------------------------------------------------
# Mode 2: feedback (rebalance recipe + redesign flagged targets)
# --------------------------------------------------------------------------

def generate_pool_recipe(feedback, output_prefix, base_vol_ul=10.0):
    rows = [{
        'Target': r['target_id'], 'Current_Rel_Mean': f"{r['rel_mean']:.2f}x",
        'Adjustment_Factor': f"{r.get('recommended_adjustment', 1.0):.2f}x", 'Action': r['action'],
        'Recommended_Vol_uL': f"{base_vol_ul * r.get('recommended_adjustment', 1.0):.2f}",
    } for r in feedback.get('target_recommendations', [])]
    if not rows:
        log.info("No recommendations in feedback; no recipe written.")
        return
    path = f"{output_prefix}_rebalancing_recipe.csv"
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    log.info("Rebalancing recipe: %s", path)


def run_feedback_mode(args):
    try:
        with open(args.feedback) as f:
            feedback = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        log.error("Cannot read feedback file '%s': %s", args.feedback, e)
        sys.exit(1)

    generate_pool_recipe(feedback, args.output_prefix)
    redesign = sorted({target_gene_id(r['target_id']) for r in feedback.get('target_recommendations', [])
                       if "REDESIGN" in r.get('action', '').upper()})
    if not redesign:
        log.info("No targets flagged for redesign.")
        return
    log.info("Redesigning %d flagged target(s) with --rescue: %s", len(redesign), ', '.join(redesign))
    target_file = f"{args.output_prefix}_redesign_targets.txt"
    with open(target_file, 'w') as f:
        f.write("\n".join(redesign) + "\n")
    design_args = argparse.Namespace(**vars(args))
    design_args.target_file = target_file
    design_args.output_prefix = f"{args.output_prefix}_redesign"
    design_args.rescue = True
    run_design_mode(design_args)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description="Design and pack multiplex PCR primer panels.")
    g = p.add_argument_group('Design')
    g.add_argument('--genome', default='common_refs/ecoli_genome.fna', help="Reference FASTA.")
    g.add_argument('--gff', default='common_refs/genomic.gff', help="GFF3 annotation.")
    g.add_argument('--target-file', help="Text file, one target gene ID per line.")
    g.add_argument('--blast-db', default='ncbi_data/ecoli_db', help="BLAST database prefix (built if missing).")
    g.add_argument('--num-candidates', type=int, default=25, help="Primer pairs designed per target.")
    g.add_argument('--monte-carlo-iters', type=int, default=20, help="Packing restarts.")
    g.add_argument('--seed', type=int, default=0, help="Seed for reproducible packing.")
    g.add_argument('--force-single-pool', action='store_true',
                   help="Keep only the largest compatible set; list targets that did not fit.")
    g.add_argument('--force-sparse', action='store_true', help="Never allow overlapping amplicons in one pool.")
    g.add_argument('--force-multiprime', action='store_true', help="Allow multi-copy primers (e.g. 16S rRNA).")
    g.add_argument('--rescue', action='store_true', help="Shift the target window if no pair passes.")
    g.add_argument('--add-hairpin-clamp', action='store_true', help="Prepend a 5' hairpin clamp (experimental).")
    t = p.add_argument_group('Thresholds (kcal/mol)')
    t.add_argument('--max-self-dimer-dg', type=float, default=CFG['self_dimer'])
    t.add_argument('--max-heterodimer-dg', type=float, default=CFG['hetero_dimer'])
    t.add_argument('--max-end-stability-dg', type=float, default=CFG['end_3p'])
    f = p.add_argument_group('Feedback')
    f.add_argument('--feedback', help="balancing_feedback.json from the analysis pipeline.")
    p.add_argument('--output-prefix', default='design/output/panel')
    p.add_argument('--log-file')
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(args.log_file)
    CFG.update(self_dimer=args.max_self_dimer_dg, hetero_dimer=args.max_heterodimer_dg,
               end_3p=args.max_end_stability_dg)
    try:
        if args.feedback:
            run_feedback_mode(args)
        elif args.target_file:
            run_design_mode(args)
        else:
            parser.error("provide --target-file (design) or --feedback (rebalance/redesign).")
    except (FileNotFoundError, RuntimeError, ValueError) as e:
        log.error("%s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
