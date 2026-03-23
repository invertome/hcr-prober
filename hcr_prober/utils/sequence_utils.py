# hcr_prober/utils/sequence_utils.py
import re, random
from Bio.Seq import Seq
from loguru import logger
IUPAC_MAP = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C']}
def reverse_complement(seq_str): return str(Seq(seq_str).reverse_complement())
def has_homopolymer(seq_str, max_len=4): return re.search(f'([ACGT])\\1{{{max_len},}}', seq_str, re.IGNORECASE) is not None
def resolve_iupac_spacer(spacer_str): return ''.join([random.choice(IUPAC_MAP.get(c.upper(), [c])) for c in spacer_str])
def parse_mask_regions(region_str):
    if not region_str: return []
    regions = []
    for part in region_str.split(','):
        try:
            start, end = map(int, part.strip().split('-'))
            regions.append((start - 1, end))
        except ValueError: logger.warning(f'Could not parse masking region \'{part}\'. Ignoring.')
    return regions