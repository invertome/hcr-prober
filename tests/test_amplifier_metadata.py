"""Tests for amplifier JSON file naming and citation metadata.

Audit 1.5. The two JSON config files were named after misleading
papers: HCR_v3_Wang_2020.json actually held the Choi 2014/2018 set
(B1-B5), and HCR_v2_Choi_2018.json actually held the Wang 2020 set
(B7/B9). After the rename:

  HCR_v3_Choi_2018.json -> B1-B5
  HCR_v3_Wang_2020.json -> B7-B9 (and any future Wang amplifiers)

Each amplifier entry must also carry a `_source` field with a
human-readable citation pointing back to the SI table that defines
its split-initiator handles.
"""
import os
import json
import pytest


def _amp_dir():
    import hcr_prober
    return os.path.join(os.path.dirname(hcr_prober.__file__), 'config', 'amplifiers')


def test_choi_json_holds_b1_through_b5():
    p = os.path.join(_amp_dir(), 'HCR_v3_Choi_2018.json')
    assert os.path.exists(p), f'HCR_v3_Choi_2018.json not found at {p}'
    with open(p) as fh:
        data = json.load(fh)
    assert set(data.keys()) >= {'B1', 'B2', 'B3', 'B4', 'B5'}
    # Wang amplifiers must NOT be in the Choi file
    assert 'B7' not in data
    assert 'B9' not in data


def test_wang_json_holds_b7_through_b9():
    p = os.path.join(_amp_dir(), 'HCR_v3_Wang_2020.json')
    assert os.path.exists(p), f'HCR_v3_Wang_2020.json not found at {p}'
    with open(p) as fh:
        data = json.load(fh)
    assert set(data.keys()) >= {'B7', 'B9'}
    # Choi amplifiers must NOT be in the Wang file
    assert 'B1' not in data
    assert 'B5' not in data


@pytest.mark.parametrize('amp_id', ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B9'])
def test_each_amplifier_has_source_citation(amp_id):
    """Every amplifier must declare its provenance in a `_source` field."""
    from hcr_prober.file_io import load_amplifiers
    import hcr_prober
    pkg_path = os.path.dirname(hcr_prober.__file__)
    amps = load_amplifiers(pkg_path)
    assert amp_id in amps, f'Amplifier {amp_id} not loaded'
    src = amps[amp_id].get('_source', '')
    assert src, f'Amplifier {amp_id} missing _source citation'
    # citation should mention an author and a year
    assert any(year in src for year in ('2014', '2018', '2020')), (
        f'_source for {amp_id} should reference a publication year, got: {src}'
    )


def test_amplifiers_still_loadable_after_rename():
    """The rename + new metadata must not break load_amplifiers."""
    from hcr_prober.file_io import load_amplifiers
    import hcr_prober
    amps = load_amplifiers(os.path.dirname(hcr_prober.__file__))
    assert {'B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B9'}.issubset(set(amps.keys()))
    for amp_id, amp in amps.items():
        assert 'up' in amp and 'dn' in amp


def test_amplifier_load_order_is_deterministic():
    """load_amplifiers iterates over glob.glob output and inserts into a dict
    in that order. Filesystem glob ordering is not guaranteed across machines,
    so we sort to keep the loaded amplifier dict order stable.
    """
    from hcr_prober.file_io import load_amplifiers
    import hcr_prober
    pkg = os.path.dirname(hcr_prober.__file__)
    a = load_amplifiers(pkg)
    b = load_amplifiers(pkg)
    assert list(a.keys()) == list(b.keys()), (
        f'amplifier load order differs across calls: {list(a.keys())} vs {list(b.keys())}'
    )
    # And we expect them to be sorted lexicographically across files
    # (HCR_v3_Choi_2018.json comes before HCR_v3_Wang_2020.json alphabetically,
    # so B1-B5 should appear before B7/B9).
    keys = list(a.keys())
    b1_idx = keys.index('B1')
    b7_idx = keys.index('B7')
    assert b1_idx < b7_idx, (
        f'B1 should appear before B7 (Choi file sorts before Wang file): {keys}'
    )
