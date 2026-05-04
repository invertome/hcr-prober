"""Regression test that locks in correct split-initiator handle orientation.

Sources for canonical sequences (verified 2026-05-03):
- Choi, Beck, Pierce 2014, "Next-Generation in Situ Hybridization Chain Reaction"
  ACS Nano 8(5):4284-4294, supplement page 43 (S6 amplifier sequences).
- Choi et al. 2018, "Third-generation in situ hybridization chain reaction"
  Development 145:dev165753 — defines the v3 split-initiator + 2-nt spacer architecture.

The v3.0 split-initiator design splits each 36-nt I1 sequence into two 18-nt halves:
  - 5' half ("up") sits at the 5' end of the upstream probe
  - 3' half ("dn") sits at the 3' end of the downstream probe
With both probes hybridized adjacently on the target, I1 is reconstituted to nucleate
the H1 hairpin and trigger HCR amplification.

Layout:
  Upstream probe   = 5'-[up_handle][up_spacer][upstream_target_arm]-3'
  Downstream probe = 5'-[downstream_target_arm][dn_spacer][dn_handle]-3'

If `finalize_probes` is ever refactored and accidentally swaps which arm gets which
handle, or attaches a handle to the wrong end, this test fails loudly.
"""
import os
import pytest
from hcr_prober.prober import finalize_probes
from hcr_prober.file_io import load_amplifiers


@pytest.fixture(scope='module')
def amplifiers():
    """Load amplifier definitions from the package config directory."""
    import hcr_prober
    pkg_path = os.path.dirname(hcr_prober.__file__)
    amps = load_amplifiers(pkg_path)
    assert amps, 'No amplifiers loaded — package config dir missing?'
    return amps


CANONICAL_HANDLES = {
    # Source: Choi 2014 ACS Nano SI page 43 (S6).
    # I1 split into first 18 nt (= up) and last 18 nt (= dn).
    'B1': {'up': 'GAGGAGGGCAGCAAACGG', 'dn': 'GAAGAGTCTTCCTTTACG'},
    'B2': {'up': 'CCTCGTAAATCCTCATCA', 'dn': 'ATCATCCAGTAAACCGCC'},
    'B3': {'up': 'GTCCCTGCCTCTATATCT', 'dn': 'CCACTCAACTTTAACCCG'},
    'B4': {'up': 'CCTCAACCTACCTCCAAC', 'dn': 'TCTCACCATATTCGCTTC'},
    'B5': {'up': 'CTCACTCCCAATCTCTAT', 'dn': 'CTACCCTACAAATCCAAT'},
}


@pytest.mark.parametrize('amp_id,handles', list(CANONICAL_HANDLES.items()))
def test_split_initiator_handle_orientation(amplifiers, amp_id, handles):
    """For each verified amplifier, finalize_probes must place handles correctly.

    - probe_up_final must begin with the canonical 'up' handle (5' end of upstream arm).
    - probe_dn_final must end with the canonical 'dn' handle (3' end of downstream arm).
    - The 25-nt target arms must appear after/before the handle+spacer respectively.
    """
    if amp_id not in amplifiers:
        pytest.skip(f'Amplifier {amp_id} not loaded in current config')
    upstream_arm = 'A' * 25  # marker — easy to spot in output
    downstream_arm = 'T' * 25
    blueprint = [{
        'pair_id': 'TestGene_cand_1',
        'probe_up_target': upstream_arm,
        'probe_dn_target': downstream_arm,
        'start_pos_on_sense': 100,
        'start_pos_rev': 50,
    }]

    class _Args:
        pass
    args = _Args()

    final = finalize_probes(blueprint, amp_id, amplifiers, args)
    assert len(final) == 1
    p = final[0]

    up_canonical = handles['up']
    dn_canonical = handles['dn']
    assert p['probe_up_final'].startswith(up_canonical), (
        f'{amp_id}: probe_up_final does not start with canonical up handle.\n'
        f"  expected prefix: {up_canonical}\n  got: {p['probe_up_final']}"
    )
    assert p['probe_dn_final'].endswith(dn_canonical), (
        f'{amp_id}: probe_dn_final does not end with canonical dn handle.\n'
        f"  expected suffix: {dn_canonical}\n  got: {p['probe_dn_final']}"
    )
    assert p['probe_up_final'].endswith(upstream_arm), (
        f'{amp_id}: upstream target arm not at 3-prime end of probe_up_final.\n'
        f"  got: {p['probe_up_final']}"
    )
    assert p['probe_dn_final'].startswith(downstream_arm), (
        f'{amp_id}: downstream target arm not at 5-prime end of probe_dn_final.\n'
        f"  got: {p['probe_dn_final']}"
    )


def test_b1_full_layout_matches_choi_2014(amplifiers):
    """End-to-end check on B1: full upstream and downstream sequences match expected layout."""
    if 'B1' not in amplifiers:
        pytest.skip('B1 amplifier not loaded')
    blueprint = [{
        'pair_id': 'TestGene_cand_1',
        'probe_up_target': 'GTTCTTCTGCTTGTCGGCCATGATA',  # 25 nt (eGFP excerpt example)
        'probe_dn_target': 'TAGACGTTGTGGCTGTTGTAGTTGT',  # 25 nt
        'start_pos_on_sense': 100,
        'start_pos_rev': 50,
    }]

    class _Args:
        pass
    args = _Args()

    final = finalize_probes(blueprint, 'B1', amplifiers, args)
    assert len(final) == 1
    p = final[0]
    up_spc = amplifiers['B1'].get('upspc', '')
    dn_spc = amplifiers['B1'].get('dnspc', '')
    expected_up = f'GAGGAGGGCAGCAAACGG{up_spc}GTTCTTCTGCTTGTCGGCCATGATA'
    expected_dn = f'TAGACGTTGTGGCTGTTGTAGTTGT{dn_spc}GAAGAGTCTTCCTTTACG'
    assert p['probe_up_final'] == expected_up, (
        f"B1 upstream layout mismatch:\n  expected: {expected_up}\n  got: {p['probe_up_final']}"
    )
    assert p['probe_dn_final'] == expected_dn, (
        f"B1 downstream layout mismatch:\n  expected: {expected_dn}\n  got: {p['probe_dn_final']}"
    )
