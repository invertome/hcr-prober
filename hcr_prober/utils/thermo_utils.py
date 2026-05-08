from functools import lru_cache
from Bio.SeqUtils import MeltingTemp as mt

# Empirical formamide constant for DNA-DNA duplexes:
# Tm decreases by ~0.65 C per 1% (v/v) formamide.
# (Casey & Davidson 1977; McConaughy et al 1969)
FORMAMIDE_TM_COEFF = 0.65

# Empirical urea constant for DNA-DNA duplexes:
# Tm decreases by ~2.25 C per molar urea.
# (Hutton 1977; consistent with the formamide constant in the same
#  body of work. Independent of and additive with formamide.)
UREA_TM_COEFF = 2.25


def calculate_gc_content(s):
    if not s:
        return 0.0
    return (s.upper().count('G') + s.upper().count('C')) / len(s) * 100


@lru_cache(maxsize=10000)
def calculate_tm(s, dnac1=25, dnac2=25, Na=825, Mg=0, dNTPs=0,
                 formamide_pct=0.0, urea_M=0.0):
    """Tm of an oligo at the given ionic conditions, with formamide and
    urea corrections.

    Tm_actual = Tm_NN(salt-corrected)
                - FORMAMIDE_TM_COEFF * formamide_pct
                - UREA_TM_COEFF * urea_M

    Defaults reflect 5xSSC ionic strength (~825 mM Na+); formamide_pct=0
    and urea_M=0 leave the salt-corrected nearest-neighbor Tm unmodified.
    For formamide HCR pass formamide_pct=50; for urea HCR pass urea_M=4.
    """
    if not s:
        return 0.0
    tm = mt.Tm_NN(s, dnac1=dnac1, dnac2=dnac2, Na=Na, Mg=Mg, dNTPs=dNTPs)
    return tm - FORMAMIDE_TM_COEFF * formamide_pct - UREA_TM_COEFF * urea_M
