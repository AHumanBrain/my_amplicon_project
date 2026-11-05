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
import datetime
import hashlib
import json
import logging
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
import tomllib
from functools import lru_cache

import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

from core.biophysics import BioConfig, get_end_stability_dg, get_heterodimer_dg, get_homodimer_dg

log = logging.getLogger("design")

# Dimer thresholds (kcal/mol); overridden from the CLI and re-applied in workers.
CFG = {
    'self_dimer': BioConfig.MAX_SELF_DIMER_DG,
    'hetero_dimer': BioConfig.MAX_HETERODIMER_DG,
    'end_3p': BioConfig.END_STABILITY_DG_THRESHOLD,
    'max_cross': 1000,   # bp; facing primers of different targets in one pool must be farther apart
}
RESCUE_SHIFTS = [0, -25, 25, -45, 45]  # bp window shifts tried with --rescue
# Truncated Illumina tails: the 3' 20 nt of the TruSeq Read 1 / Read 2 sequencing-primer sites.
# A later index PCR completes the adapters. Reads therefore start at the gene-specific primer.
P5_TAIL = 'ACACGACGCTCTTCCGATCT'   # 5' of every forward (sense) primer
P7_TAIL = 'GACGTGTGCTCTTCCGATCT'   # 5' of every reverse (antisense) primer
TARGET_COLUMNS = ['name', 'reference', 'start', 'end', 'strand', 'sequence', 'junctions', 'span', 'avoid',
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

def revcomp(seq):
    return str(Seq(seq).reverse_complement())


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

        strand = g('strand') or '+'
        if strand not in ('+', '-'):
            return None, f"{name}: 'strand' must be + or - (got '{strand}')."
        n = len(seq)
        # User positions are 1-based within the target AS GIVEN (reference orientation).
        junction_pos = _int_list(g('junctions'))   # first base AFTER the join
        if any(not 2 <= p <= n for p in junction_pos):
            return None, f"{name}: 'junctions' must be between 2 and {n} (1-based, first base after the join)."
        span = _ranges(g('span'))
        avoid = _ranges(g('avoid'))
        if any(not 1 <= a <= b <= n for a, b in span + avoid):
            return None, f"{name}: 'span'/'avoid' ranges must lie within 1-{n}."
        if strand == '-':
            # Design is always done on the sense strand: reverse-complement and flip every position.
            seq = revcomp(seq)
            junction_pos = [n + 2 - p for p in junction_pos]
            span = [(n + 1 - b, n + 1 - a) for a, b in span]
            avoid = [(n + 1 - b, n + 1 - a) for a, b in avoid]
        # primer3 wants the 0-based index of the last base BEFORE a junction (first-after is p, so p - 2).
        junctions = [p - 2 for p in junction_pos]
        lo = int(g('amplicon_min') or default_range[0])
        hi = int(g('amplicon_max') or default_range[1])
        if not 40 <= lo <= hi:
            return None, f"{name}: need 40 <= amplicon_min <= amplicon_max (got {lo}, {hi})."
        return {
            'name': name, 'seq': seq, 'contig': contig, 'offset': offset, 'strand': strand,
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


def process_single_target(target, refs, dbs, max_hits=1, force_multiprime=False,
                          rescue=False, num_candidates=25, max_product=2000, tails=('', '')):
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
            window = str(ref[lo_i:hi_i]).upper()
            seq, offset = (revcomp(window) if target['strand'] == '-' else window), lo_i
        else:
            seq, offset = target['seq'], target['offset']
        if not seq:
            continue
        n_seq = len(seq)

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
                fwd_oligo, rev_oligo = tails[0] + fwd, tails[1] + rev
                size = results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                q = results[f'PRIMER_LEFT_{i}'][0]   # index on the designed (sense) sequence
                # Report coordinates on the reference as given; for '-' targets the sense sequence
                # runs right-to-left along the reference.
                start = offset + q if target['strand'] == '+' else offset + n_seq - q - size
                products = [f"{dbname}:{c}:{sz}" for k, (dbname, _) in enumerate(dbs)
                            for c, sz in predict_products(cache[fwd][k][1], cache[rev][k][1], max_product)]
                candidates.append({
                    'target_id': name, 'strategy': strategy_idx, 'rank': i,
                    'flags': ";".join(flags) or "OK",
                    'fwd': fwd, 'rev': rev, 'fwd_oligo': fwd_oligo, 'rev_oligo': rev_oligo,
                    'fwd_tm': results[f'PRIMER_LEFT_{i}_TM'], 'rev_tm': results[f'PRIMER_RIGHT_{i}_TM'],
                    'size': size, 'contig': target['contig'], 'start': start, 'end': start + size,
                    'strand': target['strand'], 'tstart': target['offset'], 'tend': target['end'],
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


def priming_footprints(c):
    """(lo, hi, strand) of the two designed primers on the reference (0-based half-open)."""
    fl, rl = len(c['fwd']), len(c['rev'])
    if c.get('strand', '+') == '+':
        return [(c['start'], c['start'] + fl, '+'), (c['end'] - rl, c['end'], '-')]
    return [(c['end'] - fl, c['end'], '-'), (c['start'], c['start'] + rl, '+')]


def cross_product_size(a, b):
    """Smallest product two DIFFERENT targets' designed primers could make by pairing across targets
    (a forward-strand primer of one facing a reverse-strand primer of the other). None if none."""
    best = None
    for x in priming_footprints(a):
        for y in priming_footprints(b):
            if x[2] == y[2]:
                continue
            plus, minus = (x, y) if x[2] == '+' else (y, x)
            if plus[0] < minus[1]:
                size = minus[1] - plus[0]
                best = size if best is None else min(best, size)
    return best


def is_compatible(cand, pool, check_overlap):
    """Can `cand` join `pool`? Checks amplicon overlap (optional), cross-target products and primer dimers."""
    for p in pool:
        if p['contig'] != cand['contig']:
            continue
        if check_overlap and cand['start'] < p['end'] and cand['end'] > p['start']:
            return False, "amplicon overlap"
        if CFG['max_cross']:
            size = cross_product_size(cand, p)
            if size is not None and size <= CFG['max_cross']:
                return False, f"cross-target product ({size} bp)"
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
# Layouts: tiled (<= 2 pools by position) vs sparse (fewest pools by compatibility)
# --------------------------------------------------------------------------

TILED_POOLS = 2


def _spans_overlap(items):
    """True if any two (contig, lo, hi) spans overlap (sweep per contig)."""
    last = {}
    for contig, lo, hi in sorted(items):
        if contig in last and lo < last[contig]:
            return True
        last[contig] = max(last.get(contig, hi), hi)
    return False


def detect_layout(target_to_cands):
    """'tiled' if the target windows (or the best amplicons) overlap, else 'sparse'."""
    firsts = [c[0] for c in target_to_cands.values()]
    if _spans_overlap([(c['contig'], c['tstart'], c['tend']) for c in firsts]) or has_amplicon_overlap(firsts):
        return 'tiled'
    return 'sparse'


def assign_tiled_pools(target_to_cands, n_pools=TILED_POOLS):
    """Colour overlapping target windows along the genome: A,C,.. -> pool 0 and B,D,.. -> pool 1.

    Returns ({target: pool_index}, [targets that cannot be coloured with n_pools]). Targets that
    overlap no one land in pool 0.
    """
    firsts = {t: c[0] for t, c in target_to_cands.items()}
    order = sorted(firsts, key=lambda t: (firsts[t]['contig'], firsts[t]['tstart'], firsts[t]['tend'], t))
    members = [[] for _ in range(n_pools)]
    assignment, unassigned = {}, []
    for t in order:
        a = firsts[t]
        for p in range(n_pools):
            if not any(m['contig'] == a['contig'] and a['tstart'] < m['tend'] and a['tend'] > m['tstart']
                       for m in members[p]):
                members[p].append(a)
                assignment[t] = p
                break
        else:
            unassigned.append(t)
    return assignment, unassigned


def tiling_gaps(pools):
    """Internal stretches covered by no amplicon INSERT once primers are trimmed.

    Primer sequence is trimmed before alignment, so two amplicons that overlap by less than two
    primer lengths leave a blind spot of (primer_fwd + primer_rev - overlap) bp. Returns
    [(contig, start, end), ...] (0-based half-open), ignoring the outer ends of the tiled region.
    """
    by_contig = {}
    for pool in pools:
        for c in pool:
            fl, rl = len(c['fwd']), len(c['rev'])
            lo, hi = (c['start'] + fl, c['end'] - rl) if c.get('strand', '+') == '+' else (c['start'] + rl, c['end'] - fl)
            if hi > lo:
                by_contig.setdefault(c['contig'], []).append((lo, hi))
    gaps = []
    for contig, spans in by_contig.items():
        spans.sort()
        reach = spans[0][1]
        for lo, hi in spans[1:]:
            if lo > reach:
                gaps.append((contig, reach, lo))
            reach = max(reach, hi)
    return gaps


def pack_fixed_pools(target_to_cands, assignment, n_pools=TILED_POOLS, iterations=20, seed=0):
    """Pool membership is fixed by `assignment`; Monte Carlo first-fit only picks each target's primer pair.

    Within a pool, amplicons may not overlap and all primers must be dimer/cross-product compatible.
    Returns (pools, unplaced) where unplaced are targets with no compatible candidate in their pool.
    """
    rng = random.Random(seed)
    pools, unplaced = [], []
    for p in range(n_pools):
        tids = sorted(t for t, q in assignment.items() if q == p)
        if not tids:
            continue
        genomic = sorted(tids, key=lambda t: (target_to_cands[t][0]['contig'], target_to_cands[t][0]['tstart'], t))
        best, best_missing = None, None
        for it in range(iterations):
            order = genomic if it == 0 else rng.sample(tids, len(tids))
            pool, missing = [], []
            for tid in order:
                for rank, cand in enumerate(target_to_cands[tid]):
                    if is_compatible(cand, pool, True)[0]:
                        pool.append(dict(cand, rank_used=rank))
                        break
                else:
                    missing.append(tid)
            if best is None or len(missing) < len(best_missing):
                best, best_missing = pool, missing
            if not best_missing:
                break
        pools.append(best)
        unplaced += best_missing
    return pools, sorted(unplaced)


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

CSV_FIELDS = ['pool', 'target_id', 'name', 'strand', 'flags', 'fwd_oligo', 'rev_oligo', 'fwd_primer', 'rev_primer',
              'fwd_tm', 'rev_tm', 'amplicon_size', 'contig', 'start', 'end', 'specificity_hits',
              'predicted_products']


def write_outputs(pools, prefix, failed):
    """Per pool: .csv, .bed, _trimming.fasta. Overall: _order.csv (pool,name,sequence) and failures."""
    os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
    with open(f"{prefix}_order.csv", 'w', newline='') as of:
        order = csv.writer(of)
        order.writerow(['pool', 'name', 'sequence', 'role'])
        for n, pool in enumerate(pools, 1):
            rows = [{
                'pool': n, 'target_id': c['target_id'], 'name': c['name'], 'strand': c['strand'], 'flags': c['flags'],
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
                    f.write(f"{c['contig']}\t{c['start']}\t{c['end']}\t{c['name']}\t0\t{c['strand']}\n")
            # Reads begin at the gene-specific primer (tails are 5' of it and are not
            # sequenced), so trimming uses the specific primers only. <target>_F / _R split orientations.
            with open(f"{prefix}_pool_{n}_trimming.fasta", 'w') as f:
                for c in pool:
                    f.write(f">{c['target_id']}_F\n{c['fwd']}\n>{c['target_id']}_R\n{c['rev']}\n")
            for c in pool:
                order.writerow([n, f"{c['name']}_F", c['fwd_oligo'], 'forward (sense) primer'])
                order.writerow([n, f"{c['name']}_R", c['rev_oligo'], 'reverse (antisense) primer; the RT primer for stranded RNA'])
    if failed:
        with open(f"{prefix}_failed_targets.txt", 'w') as f:
            f.write("\n".join(sorted(set(failed))) + "\n")
    log.info("Wrote %d pool(s) to %s_pool_N.* and %s_order.csv", len(pools), prefix, prefix)


# --------------------------------------------------------------------------
# Mode 1: design
# --------------------------------------------------------------------------

def resolve_tails(args):
    """(forward tail, reverse tail) from the CLI. Default: truncated Illumina P5/P7; --no-tails disables."""
    if args.no_tails:
        return ('', '')
    tails = (args.fwd_tail.upper(), args.rev_tail.upper())
    if any(set(t) - set('ACGT') for t in tails):
        raise ValueError("--fwd-tail/--rev-tail must contain only A, C, G, T.")
    return tails


def run_design_mode(args):
    refs = load_references(args.reference)
    spec_fastas = args.specificity_fasta or args.reference
    dbs = []
    for fa in spec_fastas:
        prefix = os.path.join(args.blast_dir, os.path.splitext(os.path.basename(fa))[0])
        ensure_blast_db(fa, prefix)
        dbs.append((os.path.basename(prefix), prefix))

    failed, targets = [], []
    rows = read_targets(args.targets)
    for row in rows:
        t, err = resolve_target(row, refs, (args.amplicon_min, args.amplicon_max))
        (targets.append(t) if t else failed.append(err))
    opts = dict(dbs=dbs, max_hits=args.max_hits, force_multiprime=args.force_multiprime,
                rescue=args.rescue,
                num_candidates=args.num_candidates, max_product=args.max_product,
                tails=resolve_tails(args))
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

    layout = detect_layout(target_to_cands) if args.layout == 'auto' else args.layout
    gaps = []
    CFG['max_cross'] = (args.max_cross_product if args.max_cross_product is not None
                        else (args.amplicon_max if layout == 'tiled' else 1000))
    if layout == 'tiled':
        log.info("Layout: tiled. At most %d pools by position (A,C,... | B,D,...); primers are only chosen per pool. "
                 "Cross-target products shorter than %d bp are blocked.", TILED_POOLS, CFG['max_cross'])
        assignment, unassigned = assign_tiled_pools(target_to_cands)
        pools, unplaced = pack_fixed_pools(target_to_cands, assignment, TILED_POOLS,
                                           args.monte_carlo_iters, args.seed)
        for t in unassigned:
            failed.append(f"{t}: its window overlaps more than {TILED_POOLS - 1} other window(s) at once, "
                          f"so it cannot be put in {TILED_POOLS} pools.")
        for t in unplaced:
            failed.append(f"{t}: no primer pair is compatible with the other targets in its tiled pool "
                          f"(try --num-candidates, --rescue, or relaxed thresholds).")
        gaps = tiling_gaps(pools)
        if gaps:
            log.warning("Tiling leaves %d bp uncovered after primer trimming in %d stretch(es) (neighbouring "
                        "amplicons overlap by less than two primer lengths): %s", sum(h - l for _, l, h in gaps),
                        len(gaps), ", ".join(f"{c}:{l}-{h}" for c, l, h in gaps[:6]) + (" ..." if len(gaps) > 6 else ""))
            os.makedirs(os.path.dirname(args.output_prefix) or '.', exist_ok=True)
            with open(f"{args.output_prefix}_tiling_gaps.tsv", 'w') as f:
                f.write("contig\tstart\tend\tlength\n" + "".join(f"{c}\t{l}\t{h}\t{h - l}\n" for c, l, h in gaps))
    else:
        log.info("Layout: sparse (non-overlapping). Packing into the fewest compatible pools; "
                 "cross-target products shorter than %d bp are blocked.", CFG['max_cross'])
        pools = pack_targets_first_fit(target_to_cands, args.force_sparse, args.monte_carlo_iters, args.seed)
        limit = 1 if args.force_single_pool else args.max_pools
        if limit and len(pools) > limit:
            pools = sorted(pools, key=len, reverse=True)[:limit]
            left = sorted(set(target_to_cands) - {c['target_id'] for p in pools for c in p})
            msg = f"pool limit {limit}: {len(left)} target(s) could not be placed: {', '.join(left)}"
            log.warning(msg)
            failed.append(msg)

    ranks = [c['rank_used'] for p in pools for c in p]
    log.info("%d targets -> %d pool(s) %s; candidate rank used: mean %.1f, max %d",
             sum(len(p) for p in pools), len(pools), [len(p) for p in pools], sum(ranks) / len(ranks), max(ranks))
    write_outputs(pools, args.output_prefix, failed)
    manifest = build_manifest(args, layout, pools, failed, len(rows), gaps)
    with open(f"{args.output_prefix}_manifest.json", 'w') as f:
        json.dump(manifest, f, indent=2)
    log.info("Run manifest: %s_manifest.json", args.output_prefix)
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
# Project file (TOML) and run manifest (JSON)
# --------------------------------------------------------------------------

# section -> key -> (argparse dest, kind). kind: path | paths | str | int | bool | not (inverted bool)
PROJECT_SCHEMA = {
    '': {'name': ('project_name', 'str'), 'targets': ('targets', 'path'), 'reference': ('reference', 'paths'),
         'specificity_fasta': ('specificity_fasta', 'paths'), 'blast_dir': ('blast_dir', 'path'),
         'output_prefix': ('output_prefix', 'path')},
    'design': {'layout': ('layout', 'str'), 'max_pools': ('max_pools', 'int'), 'max_hits': ('max_hits', 'int'),
               'amplicon_min': ('amplicon_min', 'int'), 'amplicon_max': ('amplicon_max', 'int'),
               'max_product': ('max_product', 'int'), 'max_cross_product': ('max_cross_product', 'int'),
               'num_candidates': ('num_candidates', 'int'), 'monte_carlo_iters': ('monte_carlo_iters', 'int'),
               'seed': ('seed', 'int'), 'rescue': ('rescue', 'bool'), 'force_sparse': ('force_sparse', 'bool'),
               'force_single_pool': ('force_single_pool', 'bool'), 'force_multiprime': ('force_multiprime', 'bool')},
    'tails': {'forward': ('fwd_tail', 'str'), 'reverse': ('rev_tail', 'str'), 'enabled': ('no_tails', 'not')},
    'thresholds': {'self_dimer_dg': ('max_self_dimer_dg', 'float'), 'heterodimer_dg': ('max_heterodimer_dg', 'float'),
                   'end_stability_dg': ('max_end_stability_dg', 'float')},
}


def load_project(path):
    """Read a project TOML into argparse defaults. Paths are relative to the project file."""
    base = os.path.dirname(os.path.abspath(path))
    with open(path, 'rb') as f:
        data = tomllib.load(f)
    defaults = {}
    for section, body in [('', {k: v for k, v in data.items() if not isinstance(v, dict)})] + \
                         [(k, v) for k, v in data.items() if isinstance(v, dict)]:
        schema = PROJECT_SCHEMA.get(section)
        if schema is None:
            raise ValueError(f"{path}: unknown section [{section}]; allowed: {[k for k in PROJECT_SCHEMA if k]}")
        for key, value in body.items():
            if key not in schema:
                raise ValueError(f"{path}: unknown key '{key}'{' in [' + section + ']' if section else ''}; "
                                 f"allowed: {sorted(schema)}")
            dest, kind = schema[key]
            ok = {'str': isinstance(value, str), 'path': isinstance(value, str), 'int': type(value) is int,
                  'float': isinstance(value, (int, float)) and not isinstance(value, bool),
                  'bool': isinstance(value, bool), 'not': isinstance(value, bool),
                  'paths': isinstance(value, list) and all(isinstance(v, str) for v in value)}[kind]
            if not ok:
                raise ValueError(f"{path}: '{key}' should be of type {kind}, got {value!r}")
            if kind == 'path':
                value = os.path.join(base, value)
            elif kind == 'paths':
                value = [os.path.join(base, v) for v in value]
            elif kind == 'not':
                value = not value
            defaults[dest] = value
    return defaults


def _sha256(path, limit=64 * 1024 * 1024):
    if os.path.getsize(path) > limit:
        return None   # too large to hash cheaply; size is still recorded
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _file_info(path):
    return {'path': os.path.abspath(path), 'size_bytes': os.path.getsize(path), 'sha256': _sha256(path)}


def build_manifest(args, layout, pools, failed, n_requested, gaps=()):
    """Machine-readable record of one design run (inputs, options, results, outputs, next steps)."""
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True,
                                cwd=here, timeout=5).stdout.strip() or None
        dirty = bool(subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], capture_output=True,
                                    text=True, cwd=here, timeout=5).stdout.strip()) if commit else None
    except (OSError, subprocess.SubprocessError):
        commit = dirty = None
    prefix = args.output_prefix
    outputs = sorted({f for n in range(1, len(pools) + 1) for f in
                      (f"{prefix}_pool_{n}.csv", f"{prefix}_pool_{n}.bed", f"{prefix}_pool_{n}_trimming.fasta")}
                     | {f"{prefix}_order.csv"} | ({f"{prefix}_failed_targets.txt"} if failed else set())
                     | ({f"{prefix}_tiling_gaps.tsv"} if gaps else set()))
    ranks = [c['rank_used'] for p in pools for c in p]
    return {
        'tool': 'design.py', 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'git_commit': commit, 'git_dirty': dirty, 'project': getattr(args, 'project', None), 'name': getattr(args, 'project_name', None),
        'layout': layout, 'max_cross_product_bp': CFG['max_cross'],
        'options': {k: getattr(args, k) for k in (
            'layout', 'max_pools', 'max_hits', 'amplicon_min', 'amplicon_max', 'max_product', 'max_cross_product',
            'num_candidates', 'monte_carlo_iters', 'seed', 'rescue', 'force_sparse', 'force_single_pool',
            'force_multiprime', 'fwd_tail', 'rev_tail', 'no_tails', 'max_self_dimer_dg', 'max_heterodimer_dg',
            'max_end_stability_dg')},
        'inputs': {'targets': _file_info(args.targets), 'reference': [_file_info(f) for f in args.reference],
                   'specificity_fasta': [_file_info(f) for f in (args.specificity_fasta or args.reference)]},
        'results': {'targets_requested': n_requested, 'targets_designed': sum(len(p) for p in pools),
                    'n_pools': len(pools),
                    'pools': [{'pool': n, 'n_targets': len(p), 'targets': sorted(c['target_id'] for c in p)}
                              for n, p in enumerate(pools, 1)],
                    'mean_candidate_rank_used': round(sum(ranks) / len(ranks), 2) if ranks else None,
                    'failed': sorted(set(failed)),
                    'tiling_gaps': [{'contig': c, 'start': lo, 'end': hi, 'length': hi - lo} for c, lo, hi in gaps]},
        'outputs': outputs,
        'next_steps': [f"python scripts/design_viewer.py --prefix {prefix} --reference {' '.join(args.reference)} "
                       f"--outdir {prefix}_view"],
    }


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description="Design and pack custom-panel multiplex PCR primers.")
    p.add_argument('--project', help="Project TOML (see panel.example.toml); command-line options override it.")
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
    a = p.add_argument_group('Adapter tails (5\' of each gene-specific primer)')
    a.add_argument('--fwd-tail', default=P5_TAIL, help="Forward primer tail (default: truncated Illumina P5/Read 1).")
    a.add_argument('--rev-tail', default=P7_TAIL, help="Reverse primer tail (default: truncated Illumina P7/Read 2).")
    a.add_argument('--no-tails', action='store_true', help="Order bare gene-specific primers.")
    g.add_argument('--max-cross-product', type=int, default=None,
                   help="Primers of different targets in one pool must not face each other closer than this (bp; 0 disables). "
                        "Default: 1000 for sparse layouts, the max amplicon size for tiled (A and C share a pool, ~2 amplicons apart).")
    g.add_argument('--layout', choices=['auto', 'tiled', 'sparse'], default='auto',
                   help="tiled: at most 2 pools by position (A,C | B,D); sparse: fewest compatible pools; "
                        "auto: tiled if target windows or amplicons overlap.")
    g.add_argument('--max-pools', type=int, default=0, help="Sparse layout: keep at most this many pools (0 = no limit).")
    t = p.add_argument_group('Thresholds (kcal/mol)')
    t.add_argument('--max-self-dimer-dg', type=float, default=CFG['self_dimer'])
    t.add_argument('--max-heterodimer-dg', type=float, default=CFG['hetero_dimer'])
    t.add_argument('--max-end-stability-dg', type=float, default=CFG['end_3p'])
    f = p.add_argument_group('Feedback')
    f.add_argument('--feedback', help="balancing_feedback.json from the analysis pipeline.")
    p.add_argument('--output-prefix', default='design/output/panel')
    p.add_argument('--log-file')
    return p


def parse_args(argv=None):
    """Parse the CLI; --project supplies defaults that explicit options override."""
    parser = build_parser()
    pre, _ = parser.parse_known_args(argv)
    if pre.project:
        parser.set_defaults(**load_project(pre.project))
    return parser, parser.parse_args(argv)


def main():
    try:
        parser, args = parse_args()
    except (OSError, ValueError, tomllib.TOMLDecodeError) as e:
        sys.exit(f"error: {e}")
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
