import os, pytest, pandas as pd

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
