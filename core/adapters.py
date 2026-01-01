class Adapters:
    """Standardized NGS tail sequences (v8.0 Truncated)."""
    
    # 20-mer tails for synthesis-ready oligos (< 60nt)
    FWD_P5_TRUNCATED = 'ACACGACGCTCTTCCGATCT' 
    REV_P7_TRUNCATED = 'GACGTGTGCTCTTCCGATCT'
    
    # Full-length tails for template generation (not for synthesis)
    FWD_RC_P5_TRUNCATED = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    REV_RC_P7_TRUNCATED = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

def format_primer(specific_seq, tail=""):
    """Combines a primer-specific sequence with an adapter tail."""
    return f"{tail}{specific_seq}"
