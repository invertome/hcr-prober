from functools import lru_cache
from Bio.SeqUtils import MeltingTemp as mt

# Empirical formamide constant for DNA-DNA duplexes:
# Tm decreases by ~0.65 C per 1% (v/v) formamide.
# (Casey & Davidson 1977; McConaughy et al 1969)
FORMAMIDE_TM_COEFF = 0.65


def calculate_gc_content(s):
    if not s:
        return 0.0
    return (s.upper().count('G') + s.upper().count('C')) / len(s) * 100


@lru_cache(maxsize=10000)
def calculate_tm(s, dnac1=25, dnac2=25, Na=825, Mg=0, dNTPs=0, formamide_pct=0.0):
    """Tm of an oligo at the given ionic conditions, with formamide correction.

    Defaults reflect HCR hyb buffer (5xSSC = ~825 mM Na+); formamide_pct=0
    so the unmodified output is the salt-corrected nearest-neighbor Tm.
    Pass formamide_pct=50 to model HCR conditions (5xSSC + 50% formamide).
    """
    if not s:
        return 0.0
    tm = mt.Tm_NN(s, dnac1=dnac1, dnac2=dnac2, Na=Na, Mg=Mg, dNTPs=dNTPs)
    return tm - FORMAMIDE_TM_COEFF * formamide_pct
