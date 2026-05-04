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
