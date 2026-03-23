import pytest
import argparse
from hcr_prober.prober import generate_thermo_candidates

def make_args(**kwargs):
    defaults = dict(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, min_gc=20, max_gc=80, min_tm=0, max_tm=100,
        max_homopolymer=4, max_gc_diff=50,
        mask_regions=None, mask_sequences=None,
        na_conc=50, mg_conc=0, dntp_conc=0, dna_conc=25,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)

TEST_SEQ = (
    "ATGCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATATCGATCGATCGATCGAATTCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    "TCGATCGATCGATCGATCG"
)

def test_thermo_candidates_have_cached_values():
    args = make_args()
    candidates, audit = generate_thermo_candidates(TEST_SEQ, args)
    if len(candidates) == 0:
        pytest.skip("No candidates passed filters for this test sequence")
    c = candidates[0]
    assert 'gc_dn' in c
    assert 'gc_up' in c
    assert 'tm_dn' in c
    assert 'tm_up' in c
    assert 'probe_up_target' in c
    assert 'probe_dn_target' in c
    assert isinstance(c['gc_dn'], float)
    assert isinstance(c['tm_dn'], float)
