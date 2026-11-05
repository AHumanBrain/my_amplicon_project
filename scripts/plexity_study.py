"""In-silico study: primer-dimer potential, capacity, off-target and uniformity vs plexity (tiled design).

Stages (each caches to --workdir):
  candidates  design 25 candidate pairs per tile with real primer3 + BLAST (NO dimer screening)
  dimers      pairwise 3'-dimer thermodynamics for the naive and the brute-force-screened primer sets
  report      combine everything into results.json
See the docstrings and the printed assumptions; nothing here is wet-lab validated.
"""
import argparse
import math
import multiprocessing
import os
import pickle
import sys

import numpy as np
import primer3

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
import design  # noqa: E402
from core.biophysics import BioConfig  # noqa: E402

R_KCAL = 1.987e-3                      # kcal / (mol K)
COND = dict(mv_conc=BioConfig.SALT_MONOVALENT, dv_conc=BioConfig.SALT_DIVALENT,
            dntp_conc=BioConfig.DNTP_CONC, dna_conc=BioConfig.DNA_CONC)


def end_dg(a, b, temp_c):
    """3'-end stability of a on b, kcal/mol, at temp_c (KAPA salt conditions)."""
    return primer3.calc_end_stability(a, b, temp_c=temp_c, **COND).dg / 1000.0


def pair_dg(a, b, temp_c):
    """Worst (most negative) 3'-end dG over both directions."""
    return min(end_dg(a, b, temp_c), end_dg(b, a, temp_c))


def kd_molar(dg, temp_c):
    """Dissociation constant (M) from a duplex dG at temp_c."""
    return math.exp(dg / (R_KCAL * (temp_c + 273.15)))


def occupancy(dg, temp_c, conc_molar):
    """Fraction of a primer in this duplex at equilibrium for equal primer concentrations."""
    return conc_molar / (conc_molar + kd_molar(dg, temp_c))


def dg_threshold(temp_c, conc_molar, f_crit):
    """dG at which occupancy equals f_crit (pairs with dG below this are 'risky')."""
    return R_KCAL * (temp_c + 273.15) * math.log(conc_molar * (1 - f_crit) / f_crit)


# ---- stage 1 -------------------------------------------------------------------------
_W = {}


def _init(refs, dbs):
    _W['refs'], _W['dbs'] = refs, dbs


def _design_one(target):
    # dimer filters are disabled on purpose: screening is applied (and varied) later in this study
    design.CFG.update(self_dimer=-1e9, hetero_dimer=-1e9, end_3p=-1e9)
    return design.process_single_target(target, _W['refs'], _W['dbs'], num_candidates=25)


def stage_candidates(a):
    refs = design.load_references(a.reference)
    prefix = os.path.join(a.blast_dir, os.path.splitext(os.path.basename(a.reference[0]))[0])
    design.ensure_blast_db(a.reference[0], prefix)
    dbs = [("genome", prefix)]
    targets = []
    for i in range(a.n_targets):
        s = a.start + i * a.step
        t, err = design.resolve_target({"name": f"T{i + 1:04d}", "reference": a.contig, "start": str(s),
                                        "end": str(s + a.window - 1), "amplicon_min": str(a.amp_min),
                                        "amplicon_max": str(a.amp_max)}, refs)
        assert err is None, err
        targets.append(t)
    with multiprocessing.Pool(os.cpu_count(), initializer=_init, initargs=(refs, dbs)) as pool:
        results = pool.map(_design_one, targets, chunksize=4)
    cands = {name: c for name, c, _ in results if c}
    failed = [e for _, c, e in results if not c]
    os.makedirs(a.workdir, exist_ok=True)
    with open(os.path.join(a.workdir, "candidates.pkl"), "wb") as f:
        pickle.dump({'cands': cands, 'failed': failed, 'args': vars(a)}, f)
    print(f"{len(cands)}/{a.n_targets} tiles have candidates; {len(failed)} failed; "
          f"median candidates/tile {int(np.median([len(c) for c in cands.values()]))}")


# ---- stage 2: pairwise dimer thermodynamics ------------------------------------------
T_SETUP = 37.0          # no hot-start: dimers can extend during reaction setup
BASE_CONC = 250e-9      # M per primer, the conventional "plain" multiplex concentration
F_CRIT = occupancy(-6.0, T_SETUP, BASE_CONC)   # the tool's own limit (-6 kcal/mol at 37 C, 250 nM) ~ 0.4 %
PREFILTER = -3.0        # pairs weaker than this at 37 C stay negligible at any higher temperature


def _row_chunk(args):
    primers, rows, temp = args
    out = np.zeros((len(rows), len(primers)), dtype=np.float32)
    for k, i in enumerate(rows):
        for j in range(i, len(primers)):
            out[k, j] = pair_dg(primers[i], primers[j], temp)       # j == i is the self-dimer
    return out


def pair_matrix(primers, temp, procs=None):
    """Upper-triangular matrix of worst 3'-end dG (kcal/mol, 0 = none) for all primer pairs incl. self."""
    n = len(primers)
    order = list(range(n))
    chunks = [order[i::64] for i in range(64)]          # interleave rows so chunks have similar cost
    with multiprocessing.Pool(procs or os.cpu_count()) as pool:
        parts = pool.map(_row_chunk, [(primers, rows, temp) for rows in chunks])
    m = np.zeros((n, n), dtype=np.float32)
    for rows, part in zip(chunks, parts):
        m[rows] = part
    return m


def refine_matrix(primers, base, temp):
    """dG at `temp` only for pairs that matter at 37 C (base < PREFILTER); everything else is ~0."""
    out = np.zeros_like(base)
    idx = np.argwhere(base < PREFILTER)
    for i, j in idx:
        out[i, j] = pair_dg(primers[i], primers[j], temp)
    return out


# ---- brute-force screened selection ----------------------------------------------------
def _greedy_pool(args):
    """Greedy (genome-order) choice of one primer pair per target inside ONE fixed pool.

    A candidate is accepted if (a) its own primers (self-dimers and F/R) pass, (b) none of its primers
    has a pair with a member's primer whose 3'-end dG at `temp` is below `thr`, (c) its amplicon does
    not overlap a member's, and (d) there is no cross-target product shorter than `max_cross`.
    With repair=True, a target left unplaced is retried: if exactly ONE member blocks a candidate,
    that member is moved to one of its own alternatives (conflict-directed one-swap repair), up to
    `passes` times over the unplaced targets.
    """
    items, temp, thr, max_cross, repair, passes = args    # items: [(target, [candidate dicts])] genome order
    cache, chosen_c, chosen = {}, {}, {}
    cand_of = dict(items)

    def dg(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in cache:
            cache[key] = pair_dg(key[0], key[1], temp)
        return cache[key]

    def internal_ok(c):
        f, r = c['fwd'], c['rev']
        return dg(f, f) >= thr and dg(r, r) >= thr and dg(f, r) >= thr

    def blocks(c, m):
        """True if member candidate m prevents candidate c."""
        if c['start'] < m['end'] and c['end'] > m['start']:
            return True
        size = design.cross_product_size(c, m)
        if size is not None and size <= max_cross:
            return True
        return any(dg(x, y) < thr for x in (c['fwd'], c['rev']) for y in (m['fwd'], m['rev']))

    def blockers(c, skip=None):
        return [t for t, m in chosen_c.items() if t != skip and blocks(c, m)]

    for tid, cands in items:
        chosen[tid] = None
        for rank, c in enumerate(cands):
            if internal_ok(c) and not blockers(c):
                chosen_c[tid], chosen[tid] = c, rank
                break

    if repair:
        for _ in range(passes):
            progress = False
            for tid, cands in items:
                if chosen[tid] is not None:
                    continue
                for rank, c in enumerate(cands):
                    if not internal_ok(c):
                        continue
                    bl = blockers(c)
                    if len(bl) != 1:
                        continue
                    mt = bl[0]
                    for alt_rank, alt in enumerate(cand_of[mt]):
                        if alt is chosen_c[mt] or not internal_ok(alt) or blocks(alt, c) or blockers(alt, skip=mt):
                            continue
                        chosen_c[mt], chosen[mt] = alt, alt_rank
                        chosen_c[tid], chosen[tid] = c, rank
                        progress = True
                        break
                    if chosen[tid] is not None:
                        break
            if not progress:
                break
    return chosen


def greedy_select(cands, order, temp, conc, max_cross, repair=False, passes=3):
    """Run _greedy_pool for the two alternating tiled pools in parallel. Returns {target: rank or None}."""
    thr = dg_threshold(temp, conc, F_CRIT)
    # pool = parity of the TILE number (not the list position), so a failed tile cannot shift the alternation
    pools = [[(t, cands[t]) for t in order if (int(t[1:]) - 1) % 2 == p] for p in (0, 1)]
    with multiprocessing.Pool(2) as pool:
        parts = pool.map(_greedy_pool, [(items, temp, thr, max_cross, repair, passes) for items in pools])
    return {**parts[0], **parts[1]}


def stage_select(a):
    out_name = "selection_repair.pkl" if a.repair else "selection.pkl"
    with open(os.path.join(a.workdir, "candidates.pkl"), "rb") as f:
        cands = pickle.load(f)['cands']
    order = sorted(cands)
    max_cross = a.amp_max
    out = {'order': order}
    for name, temp, conc in (("screened_plain", T_SETUP, BASE_CONC), ("screened_rescue", a.rescue_temp, a.rescue_conc)):
        sel = greedy_select(cands, order, temp, conc, max_cross, repair=a.repair)
        placed = sum(r is not None for r in sel.values())
        print(f"{'repair' if a.repair else 'firstfit'} {name}: threshold {dg_threshold(temp, conc, F_CRIT):.2f} kcal/mol at {temp} C, {conc * 1e9:.0f} nM -> "
              f"{placed}/{len(order)} targets placed (first 8 unplaced: {[t for t in order if sel[t] is None][:8]})")
        out[name] = sel
    with open(os.path.join(a.workdir, out_name), "wb") as f:
        pickle.dump(out, f)


# ---- metrics ------------------------------------------------------------------------------
# Chemistry ladder. dimer_T = temperature at which extendable dimers matter (37 C when there is no
# hot-start, because dimers can extend during setup); anneal_T = cycling annealing temperature.
LADDER = [
    dict(key="L0", label="plain (no hot-start, 55 C anneal, 250 nM)", dimer_T=37.0, anneal_T=55.0, conc=250e-9),
    dict(key="L1", label="+ hot-start", dimer_T=55.0, anneal_T=55.0, conc=250e-9),
    dict(key="L2", label="+ stringent annealing (60 C)", dimer_T=60.0, anneal_T=60.0, conc=250e-9),
    dict(key="L3", label="+ lower primer conc (50 nM)", dimer_T=60.0, anneal_T=60.0, conc=50e-9),
    dict(key="L4", label="+ rebalanced primer conc", dimer_T=60.0, anneal_T=60.0, conc=50e-9, rebalance=True),
]
E_MAX, CYCLES = 0.95, 30        # assumed ceiling efficiency per cycle and cycle number (uniformity model)
CONC_BOUNDS = (10e-9, 250e-9)   # rebalancing may move each primer within this range


def sym_dg(upper):
    """Upper-triangular dG matrix -> symmetric matrix (0 = no interaction)."""
    return upper + upper.T - np.diag(np.diag(upper))


def sequestration(dg_sym, temp, conc_vec):
    """s_i = 1 - prod_j (1 - c_j/(c_j + Kd_ij)): fraction of primer i's 3' end tied up by any partner."""
    n = len(conc_vec)
    kd = np.full((n, n), np.inf)
    mask = dg_sym < 0
    kd[mask] = np.exp(dg_sym[mask] / (R_KCAL * (temp + 273.15)))
    partner = np.broadcast_to(conc_vec[None, :], (n, n))
    f = np.where(mask, partner / (partner + kd), 0.0)
    return 1.0 - np.exp(np.log1p(-f).sum(axis=1))


def risky_pairs(dg_upper, temp, conc):
    """Count of primer pairs (incl. self-dimers) whose occupancy at (temp, conc) is >= F_CRIT."""
    mask = dg_upper < 0
    occ = np.zeros_like(dg_upper, dtype=float)
    occ[mask] = conc / (conc + np.exp(dg_upper[mask] / (R_KCAL * (temp + 273.15))))
    return int((occ >= F_CRIT).sum())


def target_occupancy(seq, temp, conc_vec_scalar):
    """Occupancy of a primer on its perfect target at (temp, conc), from the primer-template duplex dG."""
    rc = str(design.Seq(seq).reverse_complement())
    dg = primer3.calc_heterodimer(seq, rc, temp_c=temp, **COND).dg / 1000.0
    kd = kd_molar(dg, temp)
    return conc_vec_scalar / (conc_vec_scalar + kd), kd


def uniformity(occ_f, occ_r, s_f, s_r):
    e = E_MAX * occ_f * occ_r * (1 - s_f) * (1 - s_r)
    y = (1 + e) ** CYCLES
    m = y.mean()
    return float((y >= 0.2 * m).mean()), float(np.std(np.log2(y)))


def stage_report(a):
    import json
    import time
    cands = pickle.load(open(os.path.join(a.workdir, "candidates.pkl"), "rb"))['cands']
    sels = {tag: pickle.load(open(os.path.join(a.workdir, fn), "rb"))
            for tag, fn in (("firstfit", "selection.pkl"), ("repair", "selection_repair.pkl"))
            if os.path.exists(os.path.join(a.workdir, fn))}
    order = next(iter(sels.values()))['order']
    scenarios = {'naive': {t: 0 for t in order}}                      # best primer3 pair, no screening
    for tag, sel in sels.items():
        for rule in ("plain", "rescue"):
            scenarios[f"{tag}_{rule}"] = {t: r for t, r in sel[f"screened_{rule}"].items() if r is not None}
    plex = [n for n in a.plex if n <= len(order)] + ([len(order)] if len(order) not in a.plex else [])
    res = {'order_n': len(order), 'plex': plex, 'ladder': LADDER, 'F_CRIT': F_CRIT, 'E_MAX': E_MAX,
           'CYCLES': CYCLES, 'scenarios': {}}
    for name, chosen in scenarios.items():
        cache_file = os.path.join(a.workdir, f"scenario_{name}.json")
        if os.path.exists(cache_file) and not a.force:
            res['scenarios'][name] = json.load(open(cache_file))
            print(f"{name}: loaded from cache", flush=True)
            continue
        t0 = time.time()
        rows = {0: [], 1: []}                       # per pool: [(rank_in_order, tid, F, R)]
        for k, t in enumerate(order):
            if t in chosen:
                c = cands[t][chosen[t]]
                rows[(int(t[1:]) - 1) % 2].append((k, t, c['fwd'], c['rev']))
        pools = {}
        for p, lst in rows.items():
            seqs = [x for _, _, f, r in lst for x in (f, r)]
            m37 = pair_matrix(seqs, 37.0)
            pools[p] = dict(lst=lst, seqs=seqs, m={37: m37, 55: refine_matrix(seqs, m37, 55.0),
                                                   60: refine_matrix(seqs, m37, 60.0)})
        print(f"{name}: matrices done in {time.time() - t0:.0f}s ({sum(len(v['lst']) for v in pools.values())} targets)", flush=True)
        # target occupancy / Kd per primer at each annealing temperature
        kd_t = {T: {p: np.array([target_occupancy(x, T, 1.0)[1] for x in v['seqs']]) for p, v in pools.items()}
                for T in (55.0, 60.0)}
        out = []
        for n in plex:
            cutoff = order[n - 1]
            entry = {'N': n, 'placed': sum(1 for t in chosen if t <= cutoff), 'levels': {}}
            for lv in LADDER:
                seq_all, occ_f, occ_r, risky = [], [], [], 0
                sf, sr = [], []
                for p, v in pools.items():
                    keep = [i for i, (k, t, _, _) in enumerate(v['lst']) if t <= cutoff]
                    idx = np.array([j for i in keep for j in (2 * i, 2 * i + 1)], dtype=int)
                    if len(idx) == 0:
                        continue
                    up = v['m'][int(lv['dimer_T'])][np.ix_(idx, idx)]
                    risky += risky_pairs(up, lv['dimer_T'], lv['conc'])
                    kd = kd_t[lv['anneal_T']][p][idx]
                    if lv.get('rebalance'):
                        occ_star = 0.75
                        conc = np.clip(kd * occ_star / (1 - occ_star), *CONC_BOUNDS)
                    else:
                        conc = np.full(len(idx), lv['conc'])
                    occ = conc / (conc + kd)
                    s_vec = sequestration(sym_dg(up), lv['dimer_T'], conc)
                    seq_all.append(s_vec)
                    occ_f.append(occ[0::2]); occ_r.append(occ[1::2]); sf.append(s_vec[0::2]); sr.append(s_vec[1::2])
                s_cat = np.concatenate(seq_all)
                frac, sd = uniformity(np.concatenate(occ_f), np.concatenate(occ_r), np.concatenate(sf), np.concatenate(sr))
                entry['levels'][lv['key']] = {
                    'risky_pairs': risky, 'mean_seq': float(s_cat.mean()), 'p95_seq': float(np.percentile(s_cat, 95)),
                    'frac_seq_ge_5pct': float((s_cat >= 0.05).mean()), 'uniform_frac_ge_0p2': frac, 'log2_sd': sd}
            out.append(entry)
        res['scenarios'][name] = out
        with open(cache_file, "w") as f:
            json.dump(out, f)
        print(f"{name}: done in {time.time() - t0:.0f}s", flush=True)
    with open(os.path.join(a.workdir, "results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("wrote results.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="stage", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--reference", nargs="+", default=[os.path.join(ROOT, "common_refs/ecoli_genome.fna")])
    c.add_argument("--blast-dir", default="ncbi_data")
    c.add_argument("--contig", default="NC_000913.3")
    c.add_argument("--start", type=int, default=100001)
    c.add_argument("--n-targets", type=int, default=1000)
    c.add_argument("--window", type=int, default=300)
    c.add_argument("--step", type=int, default=200)
    c.add_argument("--amp-min", type=int, default=230)
    c.add_argument("--amp-max", type=int, default=300)
    c.add_argument("--workdir", required=True)
    sl = sub.add_parser("select")
    sl.add_argument("--workdir", required=True)
    sl.add_argument("--amp-max", type=int, default=300)
    sl.add_argument("--rescue-temp", type=float, default=60.0)
    sl.add_argument("--rescue-conc", type=float, default=50e-9)
    sl.add_argument("--repair", action="store_true", help="conflict-directed one-swap repair after the greedy pass")
    r = sub.add_parser("report")
    r.add_argument("--workdir", required=True)
    r.add_argument("--plex", type=int, nargs="+", default=[10, 25, 50, 100, 200, 400, 700])
    r.add_argument("--force", action="store_true", help="ignore per-scenario caches")
    a = ap.parse_args()
    {"candidates": stage_candidates, "select": stage_select, "report": stage_report}[a.stage](a)


if __name__ == "__main__":
    main()
