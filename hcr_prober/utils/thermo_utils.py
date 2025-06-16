# hcr_prober/utils/thermo_utils.py
from Bio.SeqUtils import MeltingTemp as mt
def calculate_gc_content(s): return (s.upper().count('G') + s.upper().count('C')) / len(s) * 100 if s else 0.0
def calculate_tm(s, dnac1=25, dnac2=25, Na=50, Mg=0, dNTPs=0): return mt.Tm_NN(s, dnac1=dnac1, dnac2=dnac2, Na=Na, Mg=Mg, dNTPs=dNTPs) if s else 0.0