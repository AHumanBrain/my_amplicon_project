"""Custom-panel multiplex PCR primer design and pool packing (DNA or RNA/cDNA targets).

    python design.py --targets panel.tsv --reference ref.fa      # design + pack into pools
    python design.py --feedback coverage.json --targets panel.tsv # rebalance recipe + redesign flagged targets

The user decides what to target; this tool designs the primers, checks specificity,
screens dimers and packs primers into the fewest compatible pools. See README for the
targets-table format (coordinates or raw sequence, plus optional constraints such as
splice junctions the primers must span).
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
from functools import lru_cache

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
RESCUE_SHIFTS = [0, -25, 25, -45, 45]  # bp window shifts tried with --rescue
HAIRPIN_CLAMP_MAX_LEN = 20
TARGET_COLUMNS = ['name', 'reference', 'start', 'end', 'sequence', 'junctions', 'span', 'avoid',
                  'amplicon_min', 'amplicon_max', 'max_hits']


def configure_logging(log_file=None):
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=logging.INFO, handlers=handlers,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def sanitize_name(text):
    """BED/Picard-safe name: no whitespace or parentheses."""
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(text)).strip('_')


# --------------------------------------------------------------------------
# Targets table
# --------------------------------------------------------------------------

def _int_list(text):
    return [int(x) for x in re.split(r'[,\s]+', text.strip()) if x]


def _ranges(text):
    """'10-20,35' -> [(10, 20), (35, 35)]  (1-based, inclusive)."""
    out = []
    for part in [p for p in re.split(r'[,\s]+', text.strip()) if p]:
        a, _, b = part.partition('-')
        out.append((int(a), int(b or a)))
    return out


def read_targets(path):
    """Read the targets TSV into a list of raw row dicts (strings). '#' lines are skipped."""
    with open(path, newline='') as f:
        lines = [ln for ln in f if ln.strip() and not ln.startswith('#')]
    rows = list(csv.DictReader(lines, delimiter='\t'))
    if not rows or 'name' not in rows[0]:
        raise ValueError(f"{path}: needs a header row with at least a 'name' column.")
    unknown = set(rows[0]) - set(TARGET_COLUMNS)
    if unknown:
        raise ValueError(f"{path}: unknown column(s) {sorted(unknown)}; allowed: {TARGET_COLUMNS}")
    names = [r['name'] for r in rows]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise ValueError(f"{path}: duplicate target name(s): {sorted(dupes)}")
    return rows


def write_targets(rows, path):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def resolve_target(row, refs, default_range=(150, 250)):
    """Turn a raw row into a design-ready target dict. Returns (target, error)."""
    name = row['name']
    g = lambda k: (row.get(k) or '').strip()  # noqa: E731
    try:
        seq_text = g('sequence').upper()
        if seq_text:
            seq, contig, offset = seq_text, g('reference') or name, (int(g('start')) - 1 if g('start') else 0)
        else:
            if not (g('reference') and g('start') and g('end')):
                return None, f"{name}: give either 'sequence' or 'reference' + 'start' + 'end'."
            contig = g('reference')
            if contig not in refs:
                short = contig.split('|')[-1]
                if short not in refs:
                    return None, f"{name}: reference '{contig}' not found in the reference FASTA(s)."
                contig = short
            start, end = int(g('start')), int(g('end'))
            if not (1 <= start <= end <= len(refs[contig].seq)):
                return None, f"{name}: region {start}-{end} is outside '{contig}' (length {len(refs[contig].seq)})."
            seq, offset = str(refs[contig].seq[start - 1:end]).upper(), start - 1

        n = len(seq)
        # User positions are 1-based within the target; primer3 wants 0-based index of the last
        # base BEFORE a junction (junction column = first base AFTER the join, hence p - 2).
        junctions = [p - 2 for p in _int_list(g('junctions'))]
        if any(not 0 <= j < n - 1 for j in junctions):
            return None, f"{name}: 'junctions' must be between 2 and {n} (1-based, first base after the join)."
        span = _ranges(g('span'))
        avoid = _ranges(g('avoid'))
        if any(not 1 <= a <= b <= n for a, b in span + avoid):
            return None, f"{name}: 'span'/'avoid' ranges must lie within 1-{n}."
        lo = int(g('amplicon_min') or default_range[0])
        hi = int(g('amplicon_max') or default_range[1])
        if not 40 <= lo <= hi:
            return None, f"{name}: need 40 <= amplicon_min <= amplicon_max (got {lo}, {hi})."
        return {
            'name': name, 'seq': seq, 'contig': contig, 'offset': offset,
            'coordinate_based': not seq_text, 'end': offset + n,
            'junctions': junctions, 'span': [(a - 1, b - a + 1) for a, b in span],
            'avoid': [(a - 1, b - a + 1) for a, b in avoid],
            'amp_range': (lo, hi), 'max_hits': int(g('max_hits')) if g('max_hits') else None,
        }, None
    except ValueError as e:
        return None, f"{name}: could not parse row ({e})."


def load_references(paths):
    refs = {}
    for p in paths:
        refs.update(SeqIO.to_dict(SeqIO.parse(p, "fasta")))
    return refs


# --------------------------------------------------------------------------
# Specificity databases (BLAST+)
# --------------------------------------------------------------------------

def ensure_blast_db(fasta, prefix):
    """Build a BLAST database for `fasta` at `prefix` unless an up-to-date one exists."""
    marker = f"{prefix}.nin" if os.path.exists(f"{prefix}.nin") else f"{prefix}.nhr"
    if os.path.exists(marker) and os.path.getmtime(marker) >= os.path.getmtime(fasta):
        return
    if not shutil.which("makeblastdb"):
        raise RuntimeError("`makeblastdb` not found; install NCBI BLAST+ (see environment.yml).")
    os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
    log.info("Building BLAST database '%s' from %s", prefix, fasta)
    try:
        subprocess.run(['makeblastdb', '-in', fasta, '-dbtype', 'nucl', '-out', prefix],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"makeblastdb failed: {e.stderr}") from e


def run_blast_batch(primer_seqs, blast_db_path, threads=1):
    """Exact full-length BLAST hits for many primers in ONE blastn call.

    Returns {seq: (hit_count, [(contig, sstart, send), ...])}. On any BLAST failure the count is
    None (callers must treat that as "unknown", never as "specific").
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
        log.error("BLAST check failed: %s", e)
        return {q: (None, []) for q in seqs}
    hits = {i: [] for i in range(len(seqs))}
    for line in proc.stdout.strip().splitlines():
        p = line.split('\t')
        if len(p) >= 4:
            hits[int(p[0][1:])].append((p[1], int(p[2]), int(p[3])))
    return {q: (len(hits[i]), hits[i]) for i, q in enumerate(seqs)}


def predict_products(fwd_locs, rev_locs, max_product=2000):
    """In-silico PCR from exact primer hits in one database.

    A forward-strand hit and a reverse-strand hit that face each other on the same contig give a
    product. Returns sorted [(contig, size), ...] for sizes <= max_product.
    """
    def norm(loc):
        contig, s, e = loc
        return contig, min(s, e), max(s, e), s <= e  # strand: True = plus

    products = set()
    for a in map(norm, fwd_locs):
        for b in map(norm, rev_locs):
            if a[0] != b[0] or a[3] == b[3]:
                continue
            plus, minus = (a, b) if a[3] else (b, a)
            if plus[1] < minus[2]:
                size = minus[2] - plus[1] + 1
                if size <= max_product:
                    products.add((a[0], size))
    return sorted(products)


# --------------------------------------------------------------------------
# Per-target design: primer3 -> specificity -> dimer filters
# --------------------------------------------------------------------------

def strategies(lo, hi):
    """(index, product_size_range, (tm_min, tm_max)). Strategy 0 first, the rest are fallbacks."""
    return [
        (0, (lo, hi), (BioConfig.IDEAL_TM_MIN, BioConfig.IDEAL_TM_MAX)),
        (1, (hi, hi + 100), (57.0, 63.0)),
        (2, (max(40, lo - 50), lo), (57.0, 63.0)),
        (3, (lo, hi), (55.0, 65.0)),
    ]


def design_primers_for_sequence(sequence, name, num_candidates=25, amp_range=(150, 250),
                                tm_range=(BioConfig.IDEAL_TM_MIN, BioConfig.IDEAL_TM_MAX),
                                junctions=(), span=(), avoid=()):
    seq_args = {'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': sequence}
    if junctions:
        seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = list(junctions)  # a primer must overlap one
    if span:
        seq_args['SEQUENCE_TARGET'] = [list(s) for s in span]  # the amplicon must contain one
    if avoid:
        seq_args['SEQUENCE_EXCLUDED_REGION'] = [list(s) for s in avoid]  # no primer may overlap
    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 19, 'PRIMER_MAX_SIZE': 21, 'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': tm_range[0], 'PRIMER_MAX_TM': tm_range[1],
        'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [list(amp_range)],
        'PRIMER_NUM_RETURN': num_candidates,
        # KAPA HiFi buffer: 2.5 mM Mg++, 1.2 mM total dNTP, SantaLucia 1998 salt correction
        'PRIMER_SALT_DIVALENT': BioConfig.SALT_DIVALENT, 'PRIMER_DNTP_CONC': BioConfig.DNTP_CONC,
        'PRIMER_SALT_CORRECTIONS': 1,
    }
    return primer3.design_primers(seq_args, global_args)


def passes_dimer_filters(fwd, rev):
    """Self-dimer, F/R heterodimer and 3' extension checks for one pair."""
    if min(get_homodimer_dg(fwd), get_homodimer_dg(rev)) < CFG['self_dimer']:
        return False
    if get_heterodimer_dg(fwd, rev) < CFG['hetero_dimer']:
        return False
    return min(get_end_stability_dg(fwd, rev), get_end_stability_dg(rev, fwd)) >= CFG['end_3p']


def process_single_target(target, refs, dbs, max_hits=1, force_multiprime=False, hairpin_clamp=False,
                          rescue=False, num_candidates=25, max_product=2000):
    """Design candidate pairs for one resolved target.

    `dbs` is a list of (name, blast_prefix). A primer is accepted if its exact-hit count is
    <= max_hits in EVERY database (per-target 'max_hits' overrides). Returns (name, candidates, error).
    """
    name = target['name']
    max_hits = target['max_hits'] if target['max_hits'] is not None else max_hits
    constrained = bool(target['junctions'] or target['span'] or target['avoid'])
    shifts = RESCUE_SHIFTS if (rescue and target['coordinate_based'] and not constrained) else [0]

    for shift in shifts:
        if shift:
            log.info("Rescue: shifting window %+d bp for %s", shift, name)
            ref = refs[target['contig']].seq
            lo_i, hi_i = max(0, target['offset'] + shift), min(len(ref), target['end'] + shift)
            seq, offset = str(ref[lo_i:hi_i]).upper(), lo_i
        else:
            seq, offset = target['seq'], target['offset']
        if not seq:
            continue

        candidates, cache = [], {}

        def collect(results, strategy_idx):
            n = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if n == 0:
                return
            wanted = {results[f'PRIMER_{side}_{k}_SEQUENCE'] for k in range(n) for side in ('LEFT', 'RIGHT')}
            need = wanted - cache.keys()
            if need:
                per_db = [run_blast_batch(need, prefix) for _, prefix in dbs]
                for s in need:
                    cache[s] = [d[s] for d in per_db]
            for i in range(n):
                fwd, rev = results[f'PRIMER_LEFT_{i}_SEQUENCE'], results[f'PRIMER_RIGHT_{i}_SEQUENCE']
                counts = [c for s in (fwd, rev) for c, _ in cache[s]]
                if None in counts:  # BLAST failed: unknown specificity, never accept
                    continue
                if not force_multiprime and any(c > max_hits for c in counts):
                    continue
                if not passes_dimer_filters(fwd, rev):
                    continue
                flags = [f"Rescued{shift:+d}bp"] if shift else []
                fwd_oligo, rev_oligo = fwd, rev
                if hairpin_clamp:
                    flags.append("HairpinClamp")
                    fwd_oligo = add_hairpin_clamp(fwd, fwd, max_clamp_len=HAIRPIN_CLAMP_MAX_LEN)
                    rev_oligo = add_hairpin_clamp(rev, rev, max_clamp_len=HAIRPIN_CLAMP_MAX_LEN)
                start = offset + results[f'PRIMER_LEFT_{i}'][0]
                size = results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                products = [f"{dbname}:{c}:{sz}" for k, (dbname, _) in enumerate(dbs)
                            for c, sz in predict_products(cache[fwd][k][1], cache[rev][k][1], max_product)]
                candidates.append({
                    'target_id': name, 'strategy': strategy_idx, 'rank': i,
                    'flags': ";".join(flags) or "OK",
                    'fwd': fwd, 'rev': rev, 'fwd_oligo': fwd_oligo, 'rev_oligo': rev_oligo,
                    'fwd_tm': results[f'PRIMER_LEFT_{i}_TM'], 'rev_tm': results[f'PRIMER_RIGHT_{i}_TM'],
                    'size': size, 'contig': target['contig'], 'start': start, 'end': start + size,
                    'hits': "F:" + "/".join(str(c) for c, _ in cache[fwd]) + " R:" + "/".join(str(c) for c, _ in cache[rev]),
                    'products': ";".join(products),
                    'name': f"{sanitize_name(name)}_amplicon_S{strategy_idx}_{i}",
                    'sort_key': (strategy_idx, i),
                })

        for idx, amp, tm in strategies(*target['amp_range']):
            try:
                collect(design_primers_for_sequence(
                    seq, name, num_candidates, amp, tm,
                    target['junctions'], target['span'], target['avoid']), idx)
            except (OSError, RuntimeError) as e:  # primer3 rejects unusable templates
                log.debug("primer3 strategy %d failed for %s: %s", idx, name, e)
            if candidates:
                return name, sorted(candidates, key=lambda c: c['sort_key']), None

    why = " (constraints may be too tight)" if constrained else ""
    return name, [], f"No specific primer pair found for {name}{why}."


_W = {}  # per-worker state: references and options are shipped once, not pickled per task


def _init_worker(cfg, refs, opts):
    CFG.update(cfg)
    _W['refs'], _W['opts'] = refs, opts


def _work(target):
    return process_single_target(target, _W['refs'], **_W['opts'])


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
    tids = sorted(target_to_cands)  # worker completion order is arbitrary; sort so --seed is reproducible
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
              'fwd_tm', 'rev_tm', 'amplicon_size', 'contig', 'start', 'end', 'specificity_hits',
              'predicted_products']


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
                'specificity_hits': c['hits'], 'predicted_products': c['products']} for c in pool]
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
    refs = load_references(args.reference)
    spec_fastas = args.specificity_fasta or args.reference
    dbs = []
    for fa in spec_fastas:
        prefix = os.path.join(args.blast_dir, os.path.splitext(os.path.basename(fa))[0])
        ensure_blast_db(fa, prefix)
        dbs.append((os.path.basename(prefix), prefix))

    failed, targets = [], []
    for row in read_targets(args.targets):
        t, err = resolve_target(row, refs, (args.amplicon_min, args.amplicon_max))
        (targets.append(t) if t else failed.append(err))
    opts = dict(dbs=dbs, max_hits=args.max_hits, force_multiprime=args.force_multiprime,
                hairpin_clamp=args.add_hairpin_clamp, rescue=args.rescue,
                num_candidates=args.num_candidates, max_product=args.max_product)
    with multiprocessing.Pool(os.cpu_count(), initializer=_init_worker,
                              initargs=(dict(CFG), refs, opts)) as pool:
        results = list(tqdm(pool.imap_unordered(_work, targets), total=len(targets), desc="Designing primers"))

    target_to_cands = {}
    for tid, cands, err in results:
        if cands:
            target_to_cands[tid] = cands
        elif err:
            failed.append(err)
    if not target_to_cands:
        log.error("No primers found for any target.")
        for reason in failed:
            log.warning("Failed: %s", reason)
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

def target_name_from_amplicon(amplicon_name):
    """'gyrB_amplicon_S0_3' -> 'gyrB' (sanitised form of the target name)."""
    return amplicon_name.split('_amplicon_')[0]


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
    flagged = {target_name_from_amplicon(r['target_id']) for r in feedback.get('target_recommendations', [])
               if "REDESIGN" in r.get('action', '').upper()}
    if not flagged:
        log.info("No targets flagged for redesign.")
        return
    if not args.targets:
        log.warning("%d target(s) flagged for redesign (%s) but --targets was not given; "
                    "rerun with the original targets table to redesign them.", len(flagged), ', '.join(sorted(flagged)))
        return
    rows = [r for r in read_targets(args.targets) if sanitize_name(r['name']) in flagged]
    if not rows:
        log.error("Flagged targets %s not found in %s.", sorted(flagged), args.targets)
        sys.exit(1)
    log.info("Redesigning %d flagged target(s) with --rescue: %s", len(rows), ', '.join(r['name'] for r in rows))
    redesign_file = f"{args.output_prefix}_redesign_targets.tsv"
    write_targets(rows, redesign_file)
    design_args = argparse.Namespace(**vars(args))
    design_args.targets = redesign_file
    design_args.output_prefix = f"{args.output_prefix}_redesign"
    design_args.rescue = True
    run_design_mode(design_args)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description="Design and pack custom-panel multiplex PCR primers.")
    g = p.add_argument_group('Design')
    g.add_argument('--targets', help="Targets TSV (see README): name + coordinates or sequence, optional constraints.")
    g.add_argument('--reference', nargs='+', default=['common_refs/ecoli_genome.fna'],
                   help="Reference FASTA(s) that coordinate-based targets are cut from "
                        "(genome, transcriptome, or both).")
    g.add_argument('--specificity-fasta', nargs='+',
                   help="FASTA(s) to screen primers against (default: the reference). "
                        "For RNA targets give the transcriptome AND the genome.")
    g.add_argument('--blast-dir', default='ncbi_data', help="Where BLAST databases are built.")
    g.add_argument('--max-hits', type=int, default=1,
                   help="Allowed exact full-length hits per primer in each database "
                        "(raise for isoforms; per-target 'max_hits' column overrides).")
    g.add_argument('--amplicon-min', type=int, default=150, help="Default minimum amplicon size.")
    g.add_argument('--amplicon-max', type=int, default=250, help="Default maximum amplicon size.")
    g.add_argument('--max-product', type=int, default=2000,
                   help="Report predicted off-target products up to this size.")
    g.add_argument('--num-candidates', type=int, default=25, help="Primer pairs designed per target.")
    g.add_argument('--monte-carlo-iters', type=int, default=20, help="Packing restarts.")
    g.add_argument('--seed', type=int, default=0, help="Seed for reproducible packing.")
    g.add_argument('--force-single-pool', action='store_true',
                   help="Keep only the largest compatible set; list targets that did not fit.")
    g.add_argument('--force-sparse', action='store_true', help="Never allow overlapping amplicons in one pool.")
    g.add_argument('--force-multiprime', action='store_true', help="Allow multi-copy primers (e.g. 16S rRNA).")
    g.add_argument('--rescue', action='store_true',
                   help="Shift the window of coordinate-based, unconstrained targets if no pair passes.")
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
        elif args.targets:
            run_design_mode(args)
        else:
            parser.error("provide --targets (design) or --feedback (rebalance/redesign).")
    except (FileNotFoundError, RuntimeError, ValueError) as e:
        log.error("%s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
