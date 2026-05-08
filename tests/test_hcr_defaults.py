"""Tests for HCR-aligned default values.

Defaults model 5xSSCT + 4 M urea at 37 C (urea-HCR per Choi et al. 2018,
Development 145:dev165753) — the buffer most contemporary HCR labs use,
including this maintainer's. Reported Tm and dG must reflect those
conditions or users will tune filters against numbers that don't
correspond to their experimental reality. Formamide-HCR (5xSSC + 50%
formamide) and PCR conditions remain available as alternate presets.
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
    """5xSSC(T) is approximately 825 mM Na+ (300 mM NaCl + 30 mM Na-citrate, 5x).
    The default --buffer-preset is hcr-5xssct-urea, so after preset resolution
    Na = 825."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0, (
        f'Default --na-conc should resolve to 825 mM (5xSSCT), got {args.na_conc}'
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


# --- Default preset: urea-HCR ------------------------------------------------
def test_default_preset_is_urea_hcr():
    """Default --buffer-preset is hcr-5xssct-urea (5xSSCT + 4 M urea, 37 C).
    Resolves to Na=825 mM, formamide=0, urea_M=4."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0
    assert args.formamide_pct == 0.0, (
        f'Default --formamide-pct should resolve to 0 under urea preset, '
        f'got {args.formamide_pct}'
    )
    assert args.urea_M == 4.0, (
        f'Default --urea-M should resolve to 4 (urea-HCR), got {args.urea_M}'
    )


def test_buffer_preset_pcr_yields_pcr_ionic_values():
    """--buffer-preset pcr should populate na/mg/formamide/urea as PCR
    conditions (50 mM Na, 0 Mg, 0% formamide, 0 M urea)."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'pcr'])
    apply_buffer_preset(args)
    assert args.na_conc == 50.0
    assert args.mg_conc == 0.0
    assert args.formamide_pct == 0.0
    assert args.urea_M == 0.0


def test_buffer_preset_hcr_5xssc_yields_formamide_values():
    """--buffer-preset hcr-5xssc (formamide HCR) should populate Na=825,
    formamide=50, urea=0."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'hcr-5xssc'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0
    assert args.formamide_pct == 50.0
    assert args.urea_M == 0.0


def test_buffer_preset_hcr_5xssct_urea_yields_urea_values():
    """--buffer-preset hcr-5xssct-urea (the default) should populate Na=825,
    formamide=0, urea=4."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'hcr-5xssct-urea'])
    apply_buffer_preset(args)
    assert args.na_conc == 825.0
    assert args.formamide_pct == 0.0
    assert args.urea_M == 4.0


def test_explicit_flag_overrides_buffer_preset():
    """Explicit --na-conc beats preset; other slots still come from preset."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'hcr-5xssc',
                          '--na-conc', '100'])
    apply_buffer_preset(args)
    assert args.na_conc == 100.0
    assert args.formamide_pct == 50.0
    assert args.urea_M == 0.0


def test_explicit_urea_flag_overrides_preset():
    """Explicit --urea-M wins over preset value."""
    from hcr_prober.main import apply_buffer_preset
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--buffer-preset', 'hcr-5xssct-urea',
                          '--urea-M', '8'])
    apply_buffer_preset(args)
    assert args.urea_M == 8.0
    # other slots still from preset
    assert args.na_conc == 825.0
    assert args.formamide_pct == 0.0


# --- v1.13.0: Tm filter is opt-in -------------------------------------------
def test_tm_filter_default_is_off():
    """--min-tm and --max-tm default to None so the Tm window filter is off
    unless the user opts in. Matches contemporary urea-HCR practice where
    HCR's split-initiator architecture is tolerant enough that GC and
    homopolymer constraints alone are sufficient."""
    p = _build_design_parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert args.min_tm is None, (
        f'--min-tm default should be None (filter off), got {args.min_tm}'
    )
    assert args.max_tm is None, (
        f'--max-tm default should be None (filter off), got {args.max_tm}'
    )


def test_tm_filter_skipped_when_bounds_unset():
    """generate_thermo_candidates must NOT drop probes on Tm when both
    --min-tm and --max-tm are None. A probe whose Tm sits outside
    a hypothetical 40-55 window must still pass."""
    import argparse
    from hcr_prober.prober import generate_thermo_candidates

    # 52nt window with two GC-balanced 25-nt arms whose urea-corrected Tm
    # under default urea=4 M sits at ~62-65 C — well above 55.
    seq = 'CGATCGATCGATCGATCGATCGATC' + 'AT' + 'GATCGATCGATCGATCGATCGATCG'
    args = argparse.Namespace(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, max_homopolymer=4,
        min_gc=30.0, max_gc=70.0, max_gc_diff=20.0,
        min_tm=None, max_tm=None,
        na_conc=825.0, mg_conc=0.0, dntp_conc=0.0, dna_conc=25.0,
        formamide_pct=0.0, urea_M=4.0,
        mask_regions=None, mask_sequences=None,
    )
    candidates, audit = generate_thermo_candidates(seq, args)
    assert len(candidates) >= 1, (
        f'Filter-off should keep this probe (urea-corrected Tm >55), '
        f'got {len(candidates)} candidates. Audit: {audit}'
    )


def test_tm_filter_active_when_min_tm_passed():
    """Passing --min-tm activates the filter. Same probe with min_tm=80
    (above its urea-corrected Tm) must be dropped."""
    import argparse
    from hcr_prober.prober import generate_thermo_candidates

    seq = 'CGATCGATCGATCGATCGATCGATC' + 'AT' + 'GATCGATCGATCGATCGATCGATCG'
    args = argparse.Namespace(
        window_size=52, probe_len=25, spacer_len=2,
        skip_5prime=0, max_homopolymer=4,
        min_gc=30.0, max_gc=70.0, max_gc_diff=20.0,
        min_tm=80.0, max_tm=None,
        na_conc=825.0, mg_conc=0.0, dntp_conc=0.0, dna_conc=25.0,
        formamide_pct=0.0, urea_M=4.0,
        mask_regions=None, mask_sequences=None,
    )
    candidates, audit = generate_thermo_candidates(seq, args)
    assert len(candidates) == 0, (
        f'min_tm=80 should drop this probe (Tm ~62-65 C urea-corrected), '
        f'got {len(candidates)}. Audit: {audit}'
    )


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


# --- v1.13.0: urea Tm correction --------------------------------------------
def test_urea_correction_lowers_tm_by_225_per_M():
    """Hutton 1977: urea destabilises DNA duplexes by ~2.25 C per molar.
    With 4 M urea (standard urea-HCR), Tm should drop by ~9.0 C relative
    to urea=0."""
    from hcr_prober.utils.thermo_utils import calculate_tm
    seq = 'CGATCGATCGATCGATCGATCGATC'
    tm_no_urea = calculate_tm(seq, urea_M=0.0)
    tm_with_urea = calculate_tm(seq, urea_M=4.0)
    diff = tm_no_urea - tm_with_urea
    assert abs(diff - 9.0) < 0.05, (
        f'Expected ~9.0 C drop with 4 M urea, got {diff:.3f} C'
    )


def test_formamide_and_urea_corrections_are_additive():
    """Formamide and urea corrections are independent and additive. With
    50% formamide and 4 M urea, total Tm drop = 32.5 + 9.0 = 41.5 C
    relative to no denaturant. (Edge case; real HCR uses one or the other.)"""
    from hcr_prober.utils.thermo_utils import calculate_tm
    seq = 'CGATCGATCGATCGATCGATCGATC'
    tm_neither = calculate_tm(seq, formamide_pct=0.0, urea_M=0.0)
    tm_both = calculate_tm(seq, formamide_pct=50.0, urea_M=4.0)
    diff = tm_neither - tm_both
    assert abs(diff - 41.5) < 0.1, (
        f'Expected ~41.5 C drop (32.5 formamide + 9.0 urea), got {diff:.3f} C'
    )
