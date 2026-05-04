"""Tests for --dry-run.

Phase 5.2. --dry-run lets users tune thermodynamic / GC / Tm filter
parameters without paying the cost of running BLAST. The audit funnel
still produces, so users see how many candidates would survive into
the BLAST stage.
"""
import os
import sys
import argparse
import pytest


def test_verbose_and_quiet_flags_exist():
    """Phase 5.3: --verbose and --quiet must be CLI flags on the design parser."""
    from hcr_prober.main import add_shared_design_args
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd')
    d = sub.add_parser('design')
    d.add_argument('-i', '--input')
    d.add_argument('-o', '--output-dir')
    add_shared_design_args(d)
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert args.verbose is False
    assert args.quiet is False
    args2 = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1', '--verbose'])
    assert args2.verbose is True
    args3 = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1', '--quiet'])
    assert args3.quiet is True


def test_dry_run_flag_exists_and_defaults_false():
    from hcr_prober.main import add_shared_design_args
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd')
    d = sub.add_parser('design')
    d.add_argument('-i', '--input')
    d.add_argument('-o', '--output-dir')
    add_shared_design_args(d)
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert hasattr(args, 'dry_run')
    assert args.dry_run is False
    args2 = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1', '--dry-run'])
    assert args2.dry_run is True


def test_dry_run_skips_blast_even_with_blast_ref(monkeypatch, tmp_path):
    """--dry-run + --blast-ref must NOT invoke BLAST at all. The way main.py
    achieves this is by clearing args.blast_ref when --dry-run is set,
    which prevents create_blast_db / filter_probes_by_blast from running.
    """
    blast_calls = {'count': 0}

    def fake_create_db(ref, dbp):
        blast_calls['count'] += 1
        return None

    monkeypatch.setattr('hcr_prober.blast_wrapper.create_blast_db', fake_create_db)
    monkeypatch.setattr('hcr_prober.main.check_dependencies', lambda: None)
    monkeypatch.setattr('hcr_prober.prober.finalize_probes', lambda *a, **kw: [])
    monkeypatch.setattr('hcr_prober.prober.subsample_probes', lambda probes, n: probes)
    monkeypatch.setattr('hcr_prober.file_io.write_outputs', lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.file_io.read_fasta',
                        lambda *a, **kw: {'gene1': 'ACGT' * 100})

    fake_ref = tmp_path / 'ref.fa'
    fake_ref.write_text('>x\nACGT\n')

    monkeypatch.setattr(sys, 'argv', [
        'hcr-prober', 'design',
        '-i', '/dev/null', '-o', str(tmp_path / 'out'),
        '--amplifier', 'B1',
        '--blast-ref', str(fake_ref),
        '--dry-run',
    ])

    from hcr_prober.main import main
    main()

    assert blast_calls['count'] == 0, (
        f'create_blast_db was called {blast_calls["count"]} time(s) under '
        f'--dry-run; the flag should disable all BLAST invocations.'
    )
