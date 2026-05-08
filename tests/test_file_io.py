import argparse, os, pytest, pandas as pd

def test_csv_export_columns():
    """Verify CSV export has expected columns when probes have thermo data."""
    csv_path = '/tmp/test_details.csv'
    from hcr_prober.file_io import write_details_csv
    probes = [{
        'pair_id': 'test_pair_1', 'pair_num': 1, 'start_pos_on_sense': 100,
        'probe_dn_target': 'ATGCGATCGATCGATCGATCGATCG',
        'probe_up_target': 'TACGCTAGCTAGCTAGCTAGCTAGC',
        'gc_dn': 50.0, 'gc_up': 48.0, 'tm_dn': 60.0, 'tm_up': 58.0,
        'hairpin_dg_dn': -1.5, 'hairpin_dg_up': -1.2,
        'homodimer_dg_dn': -3.0, 'homodimer_dg_up': -2.8,
        'heterodimer_dg': -4.0,
    }]
    write_details_csv(probes, csv_path)
    assert os.path.exists(csv_path)
    df = pd.read_csv(csv_path)
    assert 'pair_id' in df.columns
    assert 'gc_dn' in df.columns
    assert 'hairpin_dg_dn' in df.columns
    assert 'heterodimer_dg' in df.columns
    assert len(df) == 1
    assert df['gc_dn'].iloc[0] == 50.0
    os.remove(csv_path)

def test_csv_export_missing_values():
    """CSV should handle probes without structure data gracefully."""
    csv_path = '/tmp/test_details_missing.csv'
    from hcr_prober.file_io import write_details_csv
    probes = [{
        'pair_id': 'test_pair_1', 'pair_num': 1, 'start_pos_on_sense': 100,
        'probe_dn_target': 'ATGCGATCGATCGATCGATCGATCG',
        'probe_up_target': 'TACGCTAGCTAGCTAGCTAGCTAGC',
        'gc_dn': 50.0, 'gc_up': 48.0, 'tm_dn': 60.0, 'tm_up': 58.0,
        # No structure values
    }]
    write_details_csv(probes, csv_path)
    df = pd.read_csv(csv_path, keep_default_na=False)
    assert df['hairpin_dg_dn'].iloc[0] == 'N/A'
    os.remove(csv_path)


# --- v1.13.0: summary thermo block reports urea + handles None Tm bounds ----
def _make_args(**overrides):
    base = dict(
        na_conc=825.0, mg_conc=0.0, dntp_conc=0.0, dna_conc=25.0,
        formamide_pct=0.0, urea_M=4.0,
        min_gc=40.0, max_gc=60.0, min_tm=None, max_tm=None,
        min_probe_distance=2,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def test_thermo_params_block_includes_urea_under_default():
    """Summary thermo block must report Urea: 4 M under default urea-HCR
    so the audit trail records the actual hyb buffer."""
    from hcr_prober.file_io import _thermo_params_block
    block = _thermo_params_block(_make_args())
    assert 'Urea: 4' in block, f'Urea field missing from thermo block:\n{block}'
    assert 'Formamide: 0' in block, f'Formamide field missing:\n{block}'
    assert 'Na+ conc: 825' in block


def test_thermo_params_block_reports_tm_filter_off_when_bounds_none():
    """When --min-tm and --max-tm are both None the block must report
    'Tm filter: off' rather than 'Tm Range: None-None C'."""
    from hcr_prober.file_io import _thermo_params_block
    block = _thermo_params_block(_make_args())
    assert 'Tm filter: off' in block, (
        f'Block should say "Tm filter: off" when bounds unset:\n{block}'
    )
    assert 'None' not in block, f'Block leaks None into output:\n{block}'


def test_thermo_params_block_reports_tm_range_when_bounds_set():
    """When user opts into the Tm filter the block reports the window."""
    from hcr_prober.file_io import _thermo_params_block
    block = _thermo_params_block(_make_args(min_tm=60.0, max_tm=75.0))
    assert 'Tm Range: 60' in block and '75' in block
    assert 'off' not in block.lower()
