"""Tests for the strict swapper matching logic.

Audit C8 / Phase 2.4. The swapper used `startswith` / `endswith`, which
matches any probe whose target portion happens to begin (or end) with a
known initiator sequence. The fix uses exact prefix/suffix length:
the first len(handle)+len(spacer) characters must equal handle+spacer
exactly, AND the total probe length must equal handle + spacer + 25.

Other regression test: a probe with the WRONG total length (e.g. a target
arm that is not 25 nt) should be left untouched, not silently mis-stripped.
"""
import os
import pandas as pd
from types import SimpleNamespace


def _amps():
    import hcr_prober
    from hcr_prober.file_io import load_amplifiers
    return load_amplifiers(os.path.dirname(hcr_prober.__file__))


def test_swap_skips_probes_with_unexpected_total_length(tmp_path):
    """If a probe is not handle + spacer + 25 nt, the swapper should not
    strip and rebuild it. The pre-fix code would still strip via
    startswith/endswith and produce a malformed output.
    """
    from hcr_prober.swapper import _swap_amplifiers_on_file
    amps = _amps()
    b1_up = amps['B1']['up']
    b1_upspc = amps['B1']['upspc']
    pool = 'B1_TestGene_PP1'
    correct = f'{b1_up}{b1_upspc}{"A" * 25}'  # 45 nt
    too_short = f'{b1_up}{b1_upspc}{"A" * 10}'  # 30 nt — non-standard
    in_path = tmp_path / 'in.xlsx'
    pd.DataFrame({'Pool name': [pool, pool], 'Sequence': [correct, too_short]}).to_excel(in_path, index=False)
    out_path = tmp_path / 'out.xlsx'
    args = SimpleNamespace(new_amplifier='B2')
    _swap_amplifiers_on_file(str(in_path), str(out_path), args, amps)
    out_df = pd.read_excel(out_path)
    seqs = list(out_df['Sequence'])
    correct_swapped = f'{amps["B2"]["up"]}{amps["B2"]["upspc"]}{"A" * 25}'
    assert seqs[0] == correct_swapped, (
        f'Standard 45-nt probe should swap to B2 layout. Got {seqs[0]!r}'
    )
    assert seqs[1] == too_short, (
        f'Non-standard length probe should be left untouched. Got {seqs[1]!r}'
    )


def test_swap_skips_target_arm_starting_with_initiator_sequence(tmp_path):
    """If the target arm coincidentally starts with another amplifier's `up`
    handle, the strict swapper must not mis-strip the handle. The probe is
    a valid B1 up-probe (correct 45 nt total, B1 up handle at the 5' end);
    swapping to B2 should yield exactly: B2 handle + B2 spacer + same 25-nt
    target arm — even though that target arm itself begins with B2's handle
    sequence.
    """
    from hcr_prober.swapper import _swap_amplifiers_on_file
    amps = _amps()
    b2_up = amps['B2']['up']  # CCTCGTAAATCCTCATCA
    target_arm = b2_up[:18] + 'TGCAGCT'  # 25 nt; starts with B2 up handle
    assert len(target_arm) == 25
    b1_probe_up = f'{amps["B1"]["up"]}{amps["B1"]["upspc"]}{target_arm}'  # 45 nt
    pool = 'B1_TestGene_PP1'
    in_path = tmp_path / 'in.xlsx'
    pd.DataFrame({'Pool name': [pool], 'Sequence': [b1_probe_up]}).to_excel(in_path, index=False)
    out_path = tmp_path / 'out.xlsx'
    args = SimpleNamespace(new_amplifier='B2')
    _swap_amplifiers_on_file(str(in_path), str(out_path), args, amps)
    out_seq = pd.read_excel(out_path)['Sequence'].iloc[0]
    expected = f'{amps["B2"]["up"]}{amps["B2"]["upspc"]}{target_arm}'
    assert out_seq == expected, (
        f'Swap mis-stripped because target arm shared B2 prefix.\n'
        f'  expected: {expected}\n'
        f'  got:      {out_seq}'
    )
