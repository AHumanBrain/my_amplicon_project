"""Same design as design.py but with wider primer limits, to test why SARS-CoV-2 tiles fail.
usage: relaxed_design.py MINSIZE,MAXSIZE,MINGC,MINTM,MAXTM <design.py options>  e.g. 20,30,30,57,64"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import primer3, design
from core.biophysics import BioConfig
smin, smax, gmin, tmlo, tmhi = [float(x) for x in sys.argv.pop(1).split(",")]
def relaxed(sequence, name, num_candidates=25, amp_range=(150, 250), tm_range=(0, 0), junctions=(), span=(), avoid=()):
    seq_args = {'SEQUENCE_ID': name, 'SEQUENCE_TEMPLATE': sequence}
    g = {'PRIMER_OPT_SIZE': int((smin + smax) // 2), 'PRIMER_MIN_SIZE': int(smin), 'PRIMER_MAX_SIZE': int(smax),
         'PRIMER_OPT_TM': 61.0, 'PRIMER_MIN_TM': tmlo, 'PRIMER_MAX_TM': tmhi, 'PRIMER_MIN_GC': gmin, 'PRIMER_MAX_GC': 60.0,
         'PRIMER_PRODUCT_SIZE_RANGE': [list(amp_range)], 'PRIMER_NUM_RETURN': num_candidates,
         'PRIMER_SALT_DIVALENT': BioConfig.SALT_DIVALENT, 'PRIMER_DNTP_CONC': BioConfig.DNTP_CONC, 'PRIMER_SALT_CORRECTIONS': 1}
    return primer3.design_primers(seq_args, g)
design.design_primers_for_sequence = relaxed
design.main()
