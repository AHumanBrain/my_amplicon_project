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

    A candidate is accepted if (a) none of its primers has a pair with a member's primer whose
    3'-end dG at `temp` is below `thr`, (b) its own F/R pair and self-dimers also pass, (c) its
    amplicon does not overlap a member's, and (d) no cross-target product shorter than `max_cross`.
    """
    items, temp, thr, max_cross = args          # items: [(target, [candidate dicts])] in genome order
    cache, members, chosen = {}, [], {}

    def dg(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in cache:
            cache[key] = pair_dg(key[0], key[1], temp)
        return cache[key]

    for tid, cands in items:
        chosen[tid] = None
        for rank, c in enumerate(cands):
            f, r = c['fwd'], c['rev']
            if dg(f, f) < thr or dg(r, r) < thr or dg(f, r) < thr:
                continue
            ok = True
            for m in members:
                if c['start'] < m['end'] and c['end'] > m['start']:
                    ok = False
                    break
                size = design.cross_product_size(c, m)
                if size is not None and size <= max_cross:
                    ok = False
                    break
                if any(dg(x, y) < thr for x in (f, r) for y in (m['fwd'], m['rev'])):
                    ok = False
                    break
            if ok:
                members.append(c)
                chosen[tid] = rank
                break
    return chosen


def greedy_select(cands, order, temp, conc, max_cross):
    """Run _greedy_pool for the two alternating tiled pools in parallel. Returns {target: rank or None}."""
    thr = dg_threshold(temp, conc, F_CRIT)
    pools = [[(t, cands[t]) for i, t in enumerate(order) if i % 2 == p] for p in (0, 1)]
    with multiprocessing.Pool(2) as pool:
        parts = pool.map(_greedy_pool, [(items, temp, thr, max_cross) for items in pools])
    return {**parts[0], **parts[1]}


def stage_select(a):
    with open(os.path.join(a.workdir, "candidates.pkl"), "rb") as f:
        cands = pickle.load(f)['cands']
    order = sorted(cands)
    max_cross = a.amp_max
    out = {'order': order}
    for name, temp, conc in (("screened_plain", T_SETUP, BASE_CONC), ("screened_rescue", a.rescue_temp, a.rescue_conc)):
        sel = greedy_select(cands, order, temp, conc, max_cross)
        placed = sum(r is not None for r in sel.values())
        print(f"{name}: threshold {dg_threshold(temp, conc, F_CRIT):.2f} kcal/mol at {temp} C, {conc * 1e9:.0f} nM -> "
              f"{placed}/{len(order)} targets placed (first 8 unplaced: {[t for t in order if sel[t] is None][:8]})")
        out[name] = sel
    with open(os.path.join(a.workdir, "selection.pkl"), "wb") as f:
        pickle.dump(out, f)


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
    a = ap.parse_args()
    {"candidates": stage_candidates, "select": stage_select}[a.stage](a)


if __name__ == "__main__":
    main()
