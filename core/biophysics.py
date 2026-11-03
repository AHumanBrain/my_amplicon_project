"""Thermodynamic helpers (primer3-py wrappers) calibrated for KAPA HiFi buffer."""
import logging
from functools import lru_cache

import primer3
from Bio.Seq import Seq

log = logging.getLogger(__name__)


class BioConfig:
    """Centralized configuration for biophysical thresholds."""
    # Reaction conditions (v10.0 KAPA HiFi calibrated)
    SALT_MONOVALENT = 50.0   # mM
    SALT_DIVALENT = 2.5      # mM Mg++
    DNTP_CONC = 1.2          # mM total
    DNA_CONC = 50.0          # nM

    # Thermodynamic thresholds (kcal/mol)
    END_STABILITY_DG_THRESHOLD = -6.0
    MAX_HETERODIMER_DG = -8.0
    MAX_SELF_DIMER_DG = -9.0
    IDEAL_TM_MIN = 59.0
    IDEAL_TM_MAX = 61.0

    # Hairpin clamp constraints
    HAIRPIN_STEM_TARGET_DG = -11.9


def _conditions():
    return dict(mv_conc=BioConfig.SALT_MONOVALENT, dv_conc=BioConfig.SALT_DIVALENT,
                dntp_conc=BioConfig.DNTP_CONC, dna_conc=BioConfig.DNA_CONC)


def _dg(result):
    """primer3 returns dG in cal/mol; convert to kcal/mol."""
    return result.dg / 1000.0


# Invalid sequences (e.g. containing N or too long) raise OSError/ValueError in
# primer3. Treat those as "unusable" (very unstable => rejected by thresholds)
# rather than silently reporting 0.0 ("no dimer"), and log it.
_BAD = -999.0


@lru_cache(maxsize=1_000_000)
def get_heterodimer_dg(seq1, seq2):
    """Binding energy between two sequences (symmetric, so keys are ordered)."""
    if seq2 < seq1:
        seq1, seq2 = seq2, seq1
    try:
        return _dg(primer3.calc_heterodimer(seq1, seq2, **_conditions()))
    except (OSError, ValueError, RuntimeError) as e:
        log.warning("calc_heterodimer failed for %s/%s: %s", seq1, seq2, e)
        return _BAD


@lru_cache(maxsize=1_000_000)
def _end_stability(seq1, seq2):
    try:
        return _dg(primer3.calc_end_stability(seq1, seq2, **_conditions()))
    except (OSError, ValueError, RuntimeError) as e:
        log.warning("calc_end_stability failed for %s/%s: %s", seq1, seq2, e)
        return _BAD


def get_end_stability_dg(seq1, seq2):
    """3' end stability of seq1 against seq2 (directional)."""
    return _end_stability(seq1, seq2)


@lru_cache(maxsize=1_000_000)
def get_homodimer_dg(seq):
    """Self-binding affinity."""
    try:
        return _dg(primer3.calc_homodimer(seq, **_conditions()))
    except (OSError, ValueError, RuntimeError) as e:
        log.warning("calc_homodimer failed for %s: %s", seq, e)
        return _BAD


def add_hairpin_clamp(oligo_sequence, primer_specific_seq,
                      target_stem_dg=BioConfig.HAIRPIN_STEM_TARGET_DG,
                      max_clamp_len=20):
    """Prepend a 5' clamp whose stem dG is closest to ``target_stem_dg``.

    Returns the oligo unchanged if no clamp improves on it.
    """
    rc = str(Seq(primer_specific_seq).reverse_complement())
    best_oligo, best_stem_dg = oligo_sequence, 0.0
    for i in range(1, min(max_clamp_len, len(rc)) + 1):
        clamp = rc[:i]
        stem_dg = get_heterodimer_dg(clamp, primer_specific_seq)
        if abs(stem_dg - target_stem_dg) < abs(best_stem_dg - target_stem_dg):
            best_oligo, best_stem_dg = clamp + oligo_sequence, stem_dg
        elif stem_dg < best_stem_dg:
            break
    return best_oligo
