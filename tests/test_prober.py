import pytest
import argparse
from hcr_prober.prober import generate_thermo_candidates
from hcr_prober.prober import select_spatially_diverse_probes
from hcr_prober.prober import filter_by_structure

def make_args(**kwargs):
    defaults = dict(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, min_gc=20, max_gc=80, min_tm=0, max_tm=100,
        max_homopolymer=4, max_gc_diff=50,
        mask_regions=None, mask_sequences=None,
        na_conc=50, mg_conc=0, dntp_conc=0, dna_conc=25,
        min_probe_distance=2,
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

def test_dp_selects_maximum_nonoverlapping():
    """Verify DP selects maximum non-overlapping set."""
    args = make_args(window_size=10, min_probe_distance=2)
    probes = [
        {'start_pos_rev': 0, 'pair_id': 'a'},
        {'start_pos_rev': 5, 'pair_id': 'b'},   # overlaps with a (0+10+2=12 > 5)
        {'start_pos_rev': 12, 'pair_id': 'c'},  # ok after a (0+12=12)
        {'start_pos_rev': 24, 'pair_id': 'd'},  # ok after c (12+12=24)
        {'start_pos_rev': 30, 'pair_id': 'e'},  # overlaps with d (24+12=36 > 30)
    ]
    result = select_spatially_diverse_probes(probes, args)
    assert len(result) == 3  # a, c, d is optimal
    ids = [p['pair_id'] for p in result]
    assert ids == ['a', 'c', 'd']

def test_dp_empty_input():
    args = make_args()
    assert select_spatially_diverse_probes([], args) == []

def test_dp_single_probe():
    args = make_args(window_size=10, min_probe_distance=2)
    probes = [{'start_pos_rev': 0, 'pair_id': 'a'}]
    result = select_spatially_diverse_probes(probes, args)
    assert len(result) == 1

def test_dp_all_overlapping():
    """All probes overlap — should select exactly 1."""
    args = make_args(window_size=10, min_probe_distance=2)
    probes = [
        {'start_pos_rev': 0, 'pair_id': 'a'},
        {'start_pos_rev': 3, 'pair_id': 'b'},
        {'start_pos_rev': 6, 'pair_id': 'c'},
    ]
    result = select_spatially_diverse_probes(probes, args)
    assert len(result) == 1

def test_dp_skip_case_critical():
    """Test where skipping a probe leads to more probes overall.
    This specifically tests the dp[i] = max(dp[i-1], ...) logic."""
    args = make_args(window_size=10, min_probe_distance=0)
    # Layout: a(0-10), b(5-15), c(10-20), d(15-25)
    # Taking a,c,d = 3 probes. Taking b,d = 2. Optimal = a,c,d
    # But b comes between a and c. Must skip b to get 3.
    probes = [
        {'start_pos_rev': 0, 'pair_id': 'a'},
        {'start_pos_rev': 5, 'pair_id': 'b'},
        {'start_pos_rev': 10, 'pair_id': 'c'},
        {'start_pos_rev': 15, 'pair_id': 'd'},  # only non-overlapping after b or c
        {'start_pos_rev': 20, 'pair_id': 'e'},
    ]
    result = select_spatially_diverse_probes(probes, args)
    assert len(result) == 3  # a, c, e (or a, b, d or similar 3-probe set)


def test_structure_filter_rejects_hairpin():
    """A known hairpin-forming sequence should be rejected."""
    args = make_args(max_hairpin_dg=-3.0, max_homodimer_dg=-5.0, max_heterodimer_dg=-5.0)
    candidates = [{
        'probe_dn_target': 'GCGCGCGCGCAAAAAGCGCGCGCGC',
        'probe_up_target': 'ATATAGATCGATCGATCGATCGATCG',
        'gc_dn': 50, 'gc_up': 40, 'tm_dn': 60, 'tm_up': 55,
    }]
    result = filter_by_structure(candidates, args)
    assert len(result) == 0, f"Expected hairpin rejection but got {len(result)} probes through"

def test_structure_filter_stores_dg_values():
    """Structure filter should store dG values in candidates that pass."""
    # Very permissive thresholds — everything passes
    args = make_args(max_hairpin_dg=-100.0, max_homodimer_dg=-100.0, max_heterodimer_dg=-100.0)
    candidates = [{
        'probe_dn_target': 'ATGCGATCGATCGATCGATCGATCG',
        'probe_up_target': 'TACGCTAGCTAGCTAGCTAGCTAGC',
        'gc_dn': 50, 'gc_up': 50, 'tm_dn': 60, 'tm_up': 60,
    }]
    result = filter_by_structure(candidates, args)
    assert len(result) == 1
    c = result[0]
    assert 'hairpin_dg_dn' in c
    assert 'hairpin_dg_up' in c
    assert 'homodimer_dg_dn' in c
    assert 'homodimer_dg_up' in c
    assert 'heterodimer_dg' in c
    assert isinstance(c['hairpin_dg_dn'], float)
