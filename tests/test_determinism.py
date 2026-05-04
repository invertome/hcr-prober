"""Tests for the --seed CLI flag and run-to-run determinism."""
import os
import sys
import argparse
import pytest


def test_skip_5prime_default_is_50():
    """Audit Phase 0.9: default --skip-5prime should be 50 (down from 100)
    to maximise the number of valid probe windows on shorter transcripts.
    """
    from hcr_prober.main import add_shared_design_args
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd')
    d = sub.add_parser('design')
    d.add_argument('-i', '--input')
    d.add_argument('-o', '--output-dir')
    add_shared_design_args(d)
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert args.skip_5prime == 50, (
        f'Default --skip-5prime should be 50, got {args.skip_5prime}'
    )


def test_seed_flag_seeds_random_and_numpy(monkeypatch, tmp_path):
    """Passing --seed N must call random.seed(N) and numpy.random.seed(N) before any RNG use."""
    captured_random = []
    captured_np = []

    import random as _random
    import numpy as _np
    monkeypatch.setattr(_random, 'seed', lambda s: captured_random.append(s))
    monkeypatch.setattr(_np.random, 'seed', lambda s: captured_np.append(s))

    monkeypatch.setattr('hcr_prober.main.check_dependencies', lambda: None)
    monkeypatch.setattr('hcr_prober.main.create_probe_blueprint',
                        lambda *a, **kw: (None, None, {}))
    monkeypatch.setattr('hcr_prober.blast_wrapper.create_blast_db',
                        lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.prober.finalize_probes', lambda *a, **kw: [])
    monkeypatch.setattr('hcr_prober.prober.subsample_probes', lambda probes, n: probes)
    monkeypatch.setattr('hcr_prober.file_io.write_outputs', lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.file_io.read_fasta',
                        lambda *a, **kw: {'gene1': 'A' * 200})

    monkeypatch.setattr(sys, 'argv', [
        'hcr-prober', 'design',
        '-i', '/dev/null',
        '-o', str(tmp_path / 'out'),
        '--amplifier', 'B1',
        '--seed', '4242',
    ])

    from hcr_prober.main import main
    main()

    assert 4242 in captured_random, (
        f'random.seed was not called with --seed value 4242. Got calls: {captured_random}'
    )
    assert 4242 in captured_np, (
        f'numpy.random.seed was not called with --seed value 4242. Got calls: {captured_np}'
    )


def test_seed_default_is_zero(monkeypatch, tmp_path):
    """If --seed is not given, default must be 0 (deterministic by default)."""
    captured_random = []

    import random as _random
    monkeypatch.setattr(_random, 'seed', lambda s: captured_random.append(s))
    monkeypatch.setattr('numpy.random.seed', lambda s: None)
    monkeypatch.setattr('hcr_prober.main.check_dependencies', lambda: None)
    monkeypatch.setattr('hcr_prober.main.create_probe_blueprint',
                        lambda *a, **kw: (None, None, {}))
    monkeypatch.setattr('hcr_prober.blast_wrapper.create_blast_db',
                        lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.prober.finalize_probes', lambda *a, **kw: [])
    monkeypatch.setattr('hcr_prober.prober.subsample_probes', lambda probes, n: probes)
    monkeypatch.setattr('hcr_prober.file_io.write_outputs', lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.file_io.read_fasta',
                        lambda *a, **kw: {'gene1': 'A' * 200})

    monkeypatch.setattr(sys, 'argv', [
        'hcr-prober', 'design',
        '-i', '/dev/null',
        '-o', str(tmp_path / 'out'),
        '--amplifier', 'B1',
    ])

    from hcr_prober.main import main
    main()

    assert 0 in captured_random, (
        f'Default --seed should be 0; random.seed(0) was not called. Got: {captured_random}'
    )


# --- Task 0.4: IUPAC spacer determinism --------------------------------------
def _pkg_path():
    import hcr_prober
    return os.path.dirname(hcr_prober.__file__)


def test_load_amplifiers_resolves_iupac_spacers():
    """After load_amplifiers, no spacer in any amplifier should contain IUPAC ambiguity codes.

    The IUPAC alphabet codes (W, R, Y, S, K, M, B, D, H, V, N) must be resolved
    to concrete bases (A/C/G/T) at load time so probe sequences are reproducible
    and the swap workflow can round-trip without sequence drift.
    """
    import random
    from hcr_prober.file_io import load_amplifiers
    random.seed(0)
    amps = load_amplifiers(_pkg_path())
    iupac_chars = set('WRYSKMBDHVN')
    for amp_id, amp in amps.items():
        for key in ('upspc', 'dnspc'):
            spc = amp.get(key, '')
            for c in spc.upper():
                assert c not in iupac_chars, (
                    f"After load, {amp_id}.{key}='{spc}' still contains unresolved "
                    f"IUPAC code '{c}'. resolve_iupac_spacer should run at load time."
                )


def test_finalize_probes_uses_amplifier_dict_spacer_verbatim():
    """finalize_probes must read upspc/dnspc straight from the amp dict and NOT
    re-resolve them. We verify this by injecting an amp with literal 'WW' and
    asserting it appears verbatim in the output. If finalize_probes still calls
    resolve_iupac_spacer, 'W' (= A/T) would be replaced and the literal 'WW'
    would not survive into the output.
    """
    from hcr_prober.prober import finalize_probes
    amps = {'TEST': {'up': 'GAGGAG', 'dn': 'TTTACG', 'upspc': 'WW', 'dnspc': 'WW'}}
    blueprint = [{
        'pair_id': 'g_cand_1',
        'probe_up_target': 'A' * 25,
        'probe_dn_target': 'T' * 25,
        'start_pos_on_sense': 0,
        'start_pos_rev': 0,
    }]

    class _Args:
        pass

    final = finalize_probes(blueprint, 'TEST', amps, _Args())
    assert len(final) == 1
    p = final[0]
    assert p['probe_up_final'] == 'GAGGAG' + 'WW' + 'A' * 25, (
        f"upstream not assembled verbatim from amp dict: {p['probe_up_final']}"
    )
    assert p['probe_dn_final'] == 'T' * 25 + 'WW' + 'TTTACG', (
        f"downstream not assembled verbatim from amp dict: {p['probe_dn_final']}"
    )
