"""Chain solver: pick amplicons along a region so that every base is covered after primer trimming.

Input is a large set of candidate amplicons (from sliding windows). The solver walks left to right choosing one
amplicon at a time, alternating two pools (A, B, A, B, ...). A chain is valid when
  * consecutive amplicon INSERTS (amplicon minus its two primers) overlap by at least `overlap` bp, so no base
    is lost to primer trimming, and
  * each new amplicon is compatible with every amplicon already in its pool (`compat`: no overlap, no
    cross-target product, no primer dimers).
The search is depth-first, furthest-reaching candidate first, with backtracking and a node budget per segment. If a
region cannot be bridged, the best chain so far is kept, the gap is recorded and a new segment starts after it.
Pure python, no primer3 or BLAST: candidates are dicts with 'start', 'end' (0-based half-open), 'fwd', 'rev', 'rank'.
"""


def insert_span(c):
    """(lo, hi) of the amplicon insert: what is left once both primers are trimmed."""
    return c['start'] + len(c['fwd']), c['end'] - len(c['rev'])


def _reach_first(c):
    return (-insert_span(c)[1], c.get('rank', 0), c['start'], c['fwd'], c['rev'])   # deterministic ties


def _search(first, parity, cands, pools, region_end, overlap, compat, end_slack, branch, budget):
    """Depth-first search for a chain that starts with `first` (already in pool `parity`).

    pools = [poolA, poolB]: amplicons already placed (other segments/regions). budget = [remaining compat checks].
    Options are evaluated lazily, furthest-reaching first, so a step normally costs one successful compat check;
    more are only evaluated when the search backtracks. Returns (best_chain, finished): the chain reaching
    furthest, and whether it reached the region end.
    """
    def options(chain):
        pool_idx = (parity + len(chain)) % 2
        members = pools[pool_idx] + [c for i, c in enumerate(chain) if (parity + i) % 2 == pool_idx]
        hi_last = insert_span(chain[-1])[1]
        for c in sorted((c for c in cands if insert_span(c)[0] <= hi_last - overlap and insert_span(c)[1] > hi_last),
                        key=_reach_first):
            budget[0] -= 1
            if budget[0] < 0:
                return
            if compat(c, members):
                yield c

    best = [first]
    stack = [([first], options([first]), [0])]            # (chain, lazy option iterator, options taken)
    while stack:
        chain, it, taken = stack[-1]
        if chain[-1]['end'] >= region_end - end_slack:
            return chain, True
        nxt = next(it, None) if taken[0] < branch and budget[0] >= 0 else None
        if nxt is None:
            stack.pop()
            continue
        taken[0] += 1
        new = chain + [nxt]
        if insert_span(nxt)[1] > insert_span(best[-1])[1]:
            best = new
        stack.append((new, options(new), [0]))
    return best, False


def chain_region(cands, region_start, region_end, pools, compat, overlap=10, end_slack=25, start_slack=25,
                 branch=6, max_nodes=3_000):
    """Chain amplicons across [region_start, region_end) (0-based half-open).

    `pools` = [poolA, poolB]: already-placed candidates, shared across regions and appended to here.
    `compat(cand, members)` -> bool says whether `cand` may join a pool holding `members`.
    Returns (placed, gaps): placed amplicons in genome order, each with 'pool' (0/1); gaps = [(lo, hi)] insert
    stretches no amplicon covers (the unavoidable outer primer sites at the region ends are not gaps).
    """
    cands = sorted((c for c in cands if region_start <= c['start'] and c['end'] <= region_end), key=_reach_first)
    placed, gaps = [], []
    parity, covered_to = 0, region_start      # next pool; insert coverage reaches here
    while True:
        pool_now = pools[parity] + [c for c in placed if c['pool'] == parity]
        ok = [c for c in cands if insert_span(c)[0] >= covered_to and compat(c, pool_now)]
        if not ok:
            if covered_to < region_end - end_slack:
                gaps.append((covered_to, region_end))
            break
        if not placed:   # first amplicon: start within the slack of the region start, furthest reaching
            near = [c for c in ok if c['start'] <= region_start + start_slack]
            first = sorted(near or ok, key=(_reach_first if near else lambda c: (c['start'],) + _reach_first(c)))[0]
        else:            # restart after a dead end: smallest gap, then furthest reaching
            first = sorted(ok, key=lambda c: (insert_span(c)[0],) + _reach_first(c))[0]
        lo_first = insert_span(first)[0]
        if lo_first > covered_to and (placed or first['start'] > region_start + start_slack):
            gaps.append((covered_to, lo_first))
        search_pools = [pools[i] + [c for c in placed if c['pool'] == i] for i in (0, 1)]
        best, finished = _search(first, parity, cands, search_pools, region_end, overlap, compat,
                                 end_slack, branch, [max_nodes])      # a fresh budget per segment
        for i, c in enumerate(best):
            c['pool'] = (parity + i) % 2
        placed.extend(best)
        if finished:
            break
        covered_to = insert_span(best[-1])[1]
        parity = (parity + len(best)) % 2
    for c in placed:
        pools[c['pool']].append(c)
    return placed, gaps
