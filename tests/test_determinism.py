"""Tests for the --seed CLI flag and run-to-run determinism.

Task 0.3 covers seeding only. The full IUPAC-spacer determinism proof
(Task 0.4) lives later in this file once that task lands.
"""
import sys
import pytest


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
