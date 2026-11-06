"""Genome-wide k-mer copy counts for target regions: a fast screen for regions that cannot hold a specific primer.

One streaming pass over a FASTA (any size; ~3 min for the human genome) counts how often each k-mer of the
target regions occurs in the whole file, on both strands. No per-target BLAST, no genome in memory.
numpy only. k <= 32; windows containing a non-ACGT base are ignored.
"""
import numpy as np

_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
    _LUT[_b + 32] = _i  # lowercase (soft-masked) bases count as bases
_HASH_MULT = np.uint64(11400714819323198485)
_HASH_BITS = 27  # 128 MB bit-per-byte table; ~0.06% false positives for 40k query k-mers, then verified exactly
_CHUNK = 8_000_000


def kmer_codes(arr, k):
    """2-bit packed forward k-mer code for every window of `arr` (uint8, 0-3 valid), and a validity mask.

    Built by doubling (1,2,4,8,...) so a chunk costs O(log k) numpy passes instead of O(k).
    """
    n = len(arr) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    bad = np.concatenate([[0], np.cumsum(arr > 3, dtype=np.int64)])
    valid = (bad[k:k + n] - bad[:n]) == 0
    base = np.where(arr > 3, 0, arr).astype(np.uint64)
    pow2 = {1: base}
    p = 1
    while p * 2 <= k:
        a = pow2[p]
        pow2[p * 2] = (a[:-p] << np.uint64(2 * p)) | a[p:]
        p *= 2
    codes, covered = None, 0
    for p in sorted(pow2, reverse=True):
        if k - covered >= p:
            part = pow2[p][covered:covered + n]
            codes = part if codes is None else (codes << np.uint64(2 * p)) | part
            covered += p
    return codes[:n], valid


def revcomp_code(code, k):
    """Reverse-complement of a packed k-mer code (python ints)."""
    out = 0
    for _ in range(k):
        out = (out << 2) | (3 - (code & 3))
        code >>= 2
    return out


def encode(seq):
    return _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]


def _fasta_chunks(path, k):
    """Yield uint8 arrays of up to _CHUNK bases (+ k-1 overlap) per record; windows never span records."""
    carry = np.empty(0, dtype=np.uint8)
    buf, size = [], 0

    def flush(final):
        nonlocal carry, buf, size
        arr = np.concatenate([carry] + [_LUT[np.frombuffer(b, dtype=np.uint8)] for b in buf]) if buf else carry
        buf, size = [], 0
        carry = arr[-(k - 1):] if (len(arr) >= k and not final) else np.empty(0, dtype=np.uint8)
        return arr

    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if buf:
                    yield flush(True)
                carry = np.empty(0, dtype=np.uint8)
                continue
            line = line.rstrip(b"\r\n")
            buf.append(line)
            size += len(line)
            if size >= _CHUNK:
                yield flush(False)
        if buf:
            yield flush(True)


def genome_copy_counts(fasta, sequences, k=19):
    """For each sequence, an int array: copies in `fasta` (both strands) of the k-mer starting at each position.

    Windows with a non-ACGT base get -1. A target cut from `fasta` itself counts at least 1 (its own site).
    """
    per_seq = []
    queries = set()
    for s in sequences:
        codes, valid = kmer_codes(encode(s.upper()), k)
        per_seq.append((codes, valid))
        queries.update(int(c) for c in codes[valid])
    canon = {}                                  # code or its reverse complement -> canonical id
    for c in queries:
        rc = revcomp_code(c, k)
        canon.setdefault(min(c, rc), len(canon))
    members = {}
    for c in queries:
        members[c] = canon[min(c, revcomp_code(c, k))]
        members.setdefault(revcomp_code(c, k), members[c])
    keys = np.array(sorted(members), dtype=np.uint64)
    ids = np.array([members[int(x)] for x in keys], dtype=np.int64)
    counts = np.zeros(len(canon), dtype=np.int64)
    if len(keys):
        table = np.zeros(1 << _HASH_BITS, dtype=bool)
        table[((keys * _HASH_MULT) >> np.uint64(64 - _HASH_BITS)).astype(np.int64)] = True
        for arr in _fasta_chunks(fasta, k):
            codes, valid = kmer_codes(arr, k)
            if not len(codes):
                continue
            hit = np.flatnonzero(valid & table[((codes * _HASH_MULT) >> np.uint64(64 - _HASH_BITS)).astype(np.int64)])
            if not len(hit):
                continue
            cand = codes[hit]
            pos = np.searchsorted(keys, cand)
            pos[pos == len(keys)] = 0
            ok = keys[pos] == cand
            counts += np.bincount(ids[pos[ok]], minlength=len(canon))
    out = []
    for codes, valid in per_seq:
        r = np.full(len(codes), -1, dtype=np.int64)
        for i in np.flatnonzero(valid):
            c = int(codes[i])
            r[i] = counts[canon[min(c, revcomp_code(c, k))]]
        out.append(r)
    return out


def has_primer_pair_site(copies, k, amp_min, amp_max, max_hits=1):
    """Can a forward and a reverse primer site (k-mers with <= max_hits genome copies) bracket an amplicon
    of amp_min..amp_max bp? `copies` are the per-position counts of the region (-1 = unusable)."""
    ok = (copies >= 0) & (copies <= max_hits)
    n = len(ok)
    cum = np.concatenate([[0], np.cumsum(ok)])
    for i in np.flatnonzero(ok):
        lo, hi = i + max(amp_min - k, 1), min(i + amp_max - k, n - 1)   # reverse site j; amplicon = i .. j+k
        if lo <= hi and cum[hi + 1] - cum[lo] > 0:
            return True
    return False
