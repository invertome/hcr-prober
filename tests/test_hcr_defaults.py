"""Tests for Phase 1 HCR-aligned default values.

The pre-Phase-1 defaults reflect PCR conditions (Na=50 mM, Mg=0,
no formamide). HCR hybridisation actually runs in 5xSSC + 50% formamide
at 37 C. The reported Tm and primer3 dG values must reflect those
conditions or users will tune filters against numbers that don't
correspond to their experimental reality.
"""
import argparse
import pytest


def _build_design_parser():
    from hcr_prober.main import add_shared_design_args
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd')
    d = sub.add_parser('design')
    d.add_argument('-i', '--input')
    d.add_argument('-o', '--output-dir')
    add_shared_design_args(d)
    return p


def test_na_conc_default_is_5xssc():
    """5xSSC is approximately 825 mM Na+ (300 mM NaCl + 30 mM Na-citrate, 5x)."""
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert args.na_conc == 825.0, (
        f'Default --na-conc should be 825 mM (5xSSC), got {args.na_conc}'
    )


# --- Task 1.2: HCR ionic strength threaded into primer3 ----------------------
def test_filter_by_structure_passes_ionic_strength_to_primer3(monkeypatch):
    """primer3.calc_hairpin/homodimer/heterodimer must receive mv_conc, dv_conc,
    dntp_conc, dna_conc, temp_c from args. Otherwise the dG values reflect the
    primer3 defaults (PCR-like) rather than HCR conditions."""
    import argparse
    from hcr_prober.prober import filter_by_structure

    captured = {'hairpin': [], 'homodimer': [], 'heterodimer': []}

    class _R:
        dg = -2000.0  # -2 kcal/mol -> passes default threshold

    def fake_hairpin(seq, **kw):
        captured['hairpin'].append(kw); return _R()
    def fake_homodimer(seq, **kw):
        captured['homodimer'].append(kw); return _R()
    def fake_heterodimer(s1, s2, **kw):
        captured['heterodimer'].append(kw); return _R()

    monkeypatch.setattr('primer3.calc_hairpin', fake_hairpin)
    monkeypatch.setattr('primer3.calc_homodimer', fake_homodimer)
    monkeypatch.setattr('primer3.calc_heterodimer', fake_heterodimer)

    args = argparse.Namespace(
        max_hairpin_dg=-3.0, max_homodimer_dg=-5.0, max_heterodimer_dg=-5.0,
        na_conc=825.0, mg_conc=10.0, dntp_conc=0.5, dna_conc=25.0,
    )
    candidates = [{
        'probe_dn_target': 'CGATCGATCGATCGATCGATCGATC',
        'probe_up_target': 'GATCGATCGATCGATCGATCGATCG',
    }]
    filter_by_structure(candidates, args)

    for fn_name, calls in captured.items():
        assert calls, f'primer3.calc_{fn_name} was never called'
        for kw in calls:
            assert kw.get('mv_conc') == 825.0, (
                f'primer3.calc_{fn_name} did not receive mv_conc=825; got {kw}'
            )
            assert kw.get('dv_conc') == 10.0, (
                f'primer3.calc_{fn_name} did not receive dv_conc=10; got {kw}'
            )
            assert kw.get('dntp_conc') == 0.5, (
                f'primer3.calc_{fn_name} did not receive dntp_conc=0.5; got {kw}'
            )
            assert kw.get('dna_conc') == 25.0, (
                f'primer3.calc_{fn_name} did not receive dna_conc=25; got {kw}'
            )
            assert kw.get('temp_c') == 37.0, (
                f'primer3.calc_{fn_name} did not receive temp_c=37; got {kw}'
            )
