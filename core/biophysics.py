import primer3
from functools import lru_cache
from Bio.Seq import Seq

class BioConfig:
    """Centralized configuration for biophysical thresholds."""
    # Salt & dNTPs (v10.0 KAPA HiFi Calibrated)
    SALT_DIVALENT = 2.5
    DNTP_CONC = 1.2
    
    # Thermodynamic Thresholds (kcal/mol)
    END_STABILITY_DG_THRESHOLD = -6.0
    MAX_HETERODIMER_DG = -8.0
    MAX_SELF_DIMER_DG = -9.0
    IDEAL_TM_MIN = 59.0
    IDEAL_TM_MAX = 61.0
    
    # Hairpin Clamp Constraints
    HAIRPIN_STEM_TARGET_DG = -11.9
    
@lru_cache(maxsize=10000)
def get_heterodimer_dg(seq1, seq2):
    """Calculates binding energy between two sequences."""
    try:
        # Note: In production we'd want to expose the salt parameters into primer3.calc_*
        return primer3.calc_heterodimer(seq1, seq2).dg / 1000.0
    except:
        return 0.0

@lru_cache(maxsize=10000)
def get_end_stability_dg(seq1, seq2):
    """Calculates 3' end stability (critical for polymerase initiation)."""
    try:
        return primer3.calc_end_stability(seq1, seq2).dg / 1000.0
    except:
        return 0.0

@lru_cache(maxsize=10000)
def get_homodimer_dg(seq):
    """Calculates self-binding affinity."""
    try:
        return primer3.calc_homodimer(seq).dg / 1000.0
    except:
        return 0.0

def add_hairpin_clamp(oligo_sequence, primer_specific_seq, 
                      target_stem_dg=BioConfig.HAIRPIN_STEM_TARGET_DG, 
                      max_clamp_len=20):
    """
    Adds a biological 'clamp' to the 5' end to suppress dimer formation 
    at high temperatures.
    """
    try:
        primer_specific_seq_rc = str(Seq(primer_specific_seq).reverse_complement())
        best_oligo = oligo_sequence
        best_stem_dg = 0.0

        for i in range(1, min(max_clamp_len, len(primer_specific_seq_rc)) + 1):
            clamp_seq = primer_specific_seq_rc[:i]
            stem_dg = get_heterodimer_dg(clamp_seq, primer_specific_seq)

            if abs(stem_dg - target_stem_dg) < abs(best_stem_dg - target_stem_dg):
                best_oligo = clamp_seq + oligo_sequence
                best_stem_dg = stem_dg
            elif stem_dg < best_stem_dg: 
                break
            
        return best_oligo
    except:
        return oligo_sequence
