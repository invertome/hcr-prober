"""Round-trip determinism test for the swap subcommand.

Depends on Task 0.4 (IUPAC spacers resolved once at amplifier load) plus
Task 0.3 (--seed flag). The swap workflow strips the old amplifier handle
and spacer, then prepends/appends the new amplifier's handle and spacer.
For amplifiers with degenerate IUPAC spacers (B7/B9 use 'ww'), the
resolution must be stable across two consecutive runs so that:

    original == swap(B7 -> B2 -> swap_back B2 -> B7)

Without the Task 0.4 fix the spacer would re-randomise on every swap
call and the round-trip would drift.
"""
import sys
import pandas as pd
import pytest


def _run_main_swap(argv, monkeypatch):
    monkeypatch.setattr(sys, 'argv', argv)
    from importlib import reload
    import hcr_prober.main as m
    reload(m)
    m.main()


def test_b7_to_b2_to_b7_roundtrip_is_deterministic(tmp_path, monkeypatch):
    target_dn = 'AAAAACCCCCGGGGGTTTTTAAAAA'
    target_up = 'TTTTTGGGGGCCCCCAAAAATTTTT'
    from hcr_prober.file_io import load_amplifiers
    import hcr_prober as _hp
    import os
    import random
    random.seed(0)
    amps = load_amplifiers(os.path.dirname(_hp.__file__))

    b7_up = amps['B7']['up']
    b7_dn = amps['B7']['dn']
    b7_upspc = amps['B7']['upspc']
    b7_dnspc = amps['B7']['dnspc']

    initial_dn_seq = f'{target_dn}{b7_dnspc}{b7_dn}'
    initial_up_seq = f'{b7_up}{b7_upspc}{target_up}'
    pool_name = 'B7_TestGene_PP1'
    initial_xlsx = tmp_path / 'initial_b7.xlsx'
    pd.DataFrame({
        'Pool name': [pool_name, pool_name],
        'Sequence': [initial_dn_seq, initial_up_seq],
    }).to_excel(initial_xlsx, index=False)

    out_b2_dir = tmp_path / 'swapped_to_b2'
    _run_main_swap([
        'hcr-prober', 'swap',
        '--input-probes', str(initial_xlsx),
        '--output-dir', str(out_b2_dir),
        '--new-amplifier', 'B2',
    ], monkeypatch)

    b2_xlsx = list(out_b2_dir.rglob('*.xlsx'))
    assert len(b2_xlsx) == 1, f'expected 1 swapped xlsx, got {b2_xlsx}'

    out_b7_dir = tmp_path / 'swapped_back_to_b7'
    _run_main_swap([
        'hcr-prober', 'swap',
        '--input-probes', str(b2_xlsx[0]),
        '--output-dir', str(out_b7_dir),
        '--new-amplifier', 'B7',
    ], monkeypatch)

    final_xlsx = list(out_b7_dir.rglob('*.xlsx'))
    assert len(final_xlsx) == 1
    final_df = pd.read_excel(final_xlsx[0])

    assert list(final_df['Sequence']) == [initial_dn_seq, initial_up_seq], (
        f'B7->B2->B7 round-trip drifted.\n'
        f"  initial:    {[initial_dn_seq, initial_up_seq]}\n"
        f"  after RT:   {list(final_df['Sequence'])}"
    )
