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
    """5xSSC is approximately 825 mM Na+ (300 mM NaCl + 30 mM Na-citrate, 5x).
    The default --buffer-preset is hcr-5xssc, so after preset resolution Na = 825."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0, (
        f'Default --na-conc should resolve to 825 mM (5xSSC), got {args.na_conc}'
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


# --- Task 1.3: formamide Tm correction --------------------------------------
def test_formamide_default_is_50pct():
    """HCR hybridisation buffer is typically 5xSSC + 50% formamide; default
    --buffer-preset hcr-5xssc resolves --formamide-pct to 50."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    apply_buffer_preset(args)
    assert args.formamide_pct == 50.0, (
        f'Default --formamide-pct should resolve to 50 (HCR), got {args.formamide_pct}'
    )


def test_buffer_preset_pcr_yields_pcr_ionic_values():
    """--buffer-preset pcr should populate na/mg/formamide as PCR conditions
    (50 mM Na, 0 Mg, 0% formamide), provided the user did NOT pass those
    individual flags explicitly.
    """
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'pcr'])
    apply_buffer_preset(args)
    assert args.na_conc == 50.0
    assert args.mg_conc == 0.0
    assert args.formamide_pct == 0.0


def test_buffer_preset_hcr_5xssc_yields_hcr_ionic_values():
    """--buffer-preset hcr-5xssc (the default) should populate Na=825, Mg=0,
    formamide=50."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0
    assert args.formamide_pct == 50.0


def test_explicit_flag_overrides_buffer_preset():
    """If the user passes --na-conc explicitly, it should win over the preset."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'hcr-5xssc',
                          '--na-conc', '100'])
    apply_buffer_preset(args)
    assert args.na_conc == 100.0
    # formamide should still come from preset since user did not override it
    assert args.formamide_pct == 50.0


def test_formamide_correction_lowers_tm_by_065_per_pct():
    """The standard rule of thumb for formamide on Tm is dT = -0.65 C / 1%
    formamide. With 50% formamide, Tm should drop by ~32.5 C relative to
    formamide=0."""
    from hcr_prober.utils.thermo_utils import calculate_tm
    seq = 'CGATCGATCGATCGATCGATCGATC'
    tm_no_formamide = calculate_tm(seq, formamide_pct=0.0)
    tm_with_formamide = calculate_tm(seq, formamide_pct=50.0)
    diff = tm_no_formamide - tm_with_formamide
    assert abs(diff - 32.5) < 0.1, (
        f'Expected ~32.5 C drop with 50% formamide, got {diff:.3f} C'
    )
