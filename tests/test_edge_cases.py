"""Edge-case tests covering tricky inputs to the design and swap workflows.

Phase 3.4. Tests for:
  - Sequence shorter than the sliding-window size (52 nt by default).
  - IUPAC ambiguity codes other than N (R/Y/S/W/K/M).
  - Swap subcommand against a directory containing multiple xlsx files.
"""
import os
import shutil
import subprocess
from types import SimpleNamespace

import pandas as pd
import pytest


def test_sequence_shorter_than_window_yields_zero_candidates():
    """A sequence shorter than --window-size cannot produce any 52-mer
    windows; generate_thermo_candidates must return [] without crashing."""
    from hcr_prober.prober import generate_thermo_candidates
    args = SimpleNamespace(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, mask_regions=None, mask_sequences=None,
        max_homopolymer=4, min_gc=20.0, max_gc=80.0, max_gc_diff=20.0,
        min_tm=0.0, max_tm=100.0,
    )
    cands, audit = generate_thermo_candidates('ACGT' * 8, args)  # 32 nt < 52
    assert cands == []
    assert audit.get('initial_windows', 0) == 0


def test_iupac_ry_filtered_alongside_n():
    """The non-ACGT filter must reject R, Y, W and other IUPAC ambiguity
    codes too, not only N."""
    from hcr_prober.prober import generate_thermo_candidates
    args = SimpleNamespace(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, mask_regions=None, mask_sequences=None,
        max_homopolymer=200,  # disable
        min_gc=0.0, max_gc=100.0, max_gc_diff=200.0,
        min_tm=0.0, max_tm=200.0,
    )
    seq = 'A' * 30 + 'RYWSKM' * 5 + 'A' * 30  # IUPAC mid-region
    cands, _ = generate_thermo_candidates(seq, args)
    for c in cands:
        seq_chars = set(c['window_sequence'].upper())
        assert seq_chars <= {'A', 'C', 'G', 'T'}, (
            f'window contained non-ACGT bases: {seq_chars}'
        )


@pytest.mark.skipif(
    not (shutil.which('blastn') and shutil.which('makeblastdb')),
    reason='hcr-prober console script test',
)
def test_swap_on_directory_of_xlsx(tmp_path):
    """swap on a directory should transform every .xlsx in it."""
    from hcr_prober.file_io import load_amplifiers
    import hcr_prober
    amps = load_amplifiers(os.path.dirname(hcr_prober.__file__))
    in_dir = tmp_path / 'input_probes'
    in_dir.mkdir()

    for i, name in enumerate(['gene_a', 'gene_b']):
        b1_up = amps['B1']['up'] + amps['B1']['upspc'] + 'A' * 25
        b1_dn = 'C' * 25 + amps['B1']['dnspc'] + amps['B1']['dn']
        pool = f'B1_{name}_PP1'
        pd.DataFrame({
            'Pool name': [pool, pool],
            'Sequence': [b1_dn, b1_up],
        }).to_excel(in_dir / f'{name}_order.xlsx', index=False)

    out_dir = tmp_path / 'swapped'
    result = subprocess.run([
        'hcr-prober', 'swap',
        '--input-probes', str(in_dir),
        '--output-dir', str(out_dir),
        '--new-amplifier', 'B5',
    ], capture_output=True, text=True)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'

    swapped_files = sorted(out_dir.glob('*.xlsx'))
    assert len(swapped_files) == 2, f'expected 2 swapped files, got {len(swapped_files)}'
    for fp in swapped_files:
        df = pd.read_excel(fp)
        for seq in df['Sequence']:
            # Should now contain B5 handles
            assert (
                seq.startswith(amps['B5']['up']) or seq.endswith(amps['B5']['dn'])
            ), f'{fp.name}: B5 handle not found in {seq!r}'
