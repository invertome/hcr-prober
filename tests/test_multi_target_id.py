"""Tests for multi-target-ID support in --target-transcript-id.

Audit C8 / Phase 2.3 — caught during functional verification on real
cerebrin transcripts. The auto-derived negative reference excludes the
target ID, but a transcriptome can contain paralogs or the antisense
contig of the same gene. Without a way to mark all of those as
non-targets, the negative screen rejects every probe.

Fix: --target-transcript-id accepts one or more IDs (nargs='+'); all
of them are added to the target_ids set used by both the
'specific-id' BLAST strategy and the negative-screen exclusion list.
"""
import argparse
import pytest


def _parser():
    from hcr_prober.main import add_shared_design_args
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd')
    d = sub.add_parser('design')
    d.add_argument('-i', '--input')
    d.add_argument('-o', '--output-dir')
    add_shared_design_args(d)
    return p


def test_target_transcript_id_accepts_multiple_values():
    p = _parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1',
                          '--target-transcript-id', 'GeneA', 'GeneB', 'GeneC'])
    assert args.target_transcript_id == ['GeneA', 'GeneB', 'GeneC']


def test_target_transcript_id_default_is_empty_list():
    p = _parser()
    args = p.parse_args(['design', '-i', 'foo', '--amplifier', 'B1'])
    assert args.target_transcript_id in (None, [])


def test_negative_screen_excludes_all_target_ids(tmp_path, monkeypatch):
    """When --target-transcript-id has multiple IDs, the auto-derived negative
    reference excludes ALL of them, not just the first."""
    from hcr_prober.blast_wrapper import run_negative_screen
    from types import SimpleNamespace

    neg_ref_path_holder = {}

    def fake_run(cmd, **kw):
        try:
            oi = cmd.index('-out')
            open(cmd[oi + 1], 'w').close()
        except ValueError:
            pass
        return SimpleNamespace(returncode=0, stderr='', stdout='')

    monkeypatch.setattr('subprocess.run', fake_run)

    blast_ref = tmp_path / 'ref.fasta'
    blast_ref.write_text(
        '>GeneA\n' + 'ACGT' * 50 + '\n'
        '>GeneA_AS\n' + 'TGCA' * 50 + '\n'
        '>GeneA_paralog\n' + 'ACTG' * 50 + '\n'
        '>UnrelatedGene\n' + 'GGCC' * 50 + '\n'
    )
    args = SimpleNamespace(
        blast_ref=str(blast_ref),
        blast_negative_ref=None,
        negative_bitscore=None,
        negative_evalue=None,
        min_bitscore=75.0,
        max_evalue=1e-10,
        blast_extra_args=[],
    )
    target_ids = {'GeneA', 'GeneA_AS', 'GeneA_paralog'}

    probes = [{'pair_id': 'p1', 'probe_dn_target': 'A' * 25, 'probe_up_target': 'T' * 25}]
    run_negative_screen(probes, args, str(tmp_path), target_ids)

    neg_ref_files = list((tmp_path).glob('negative_ref.fasta'))
    assert neg_ref_files, 'auto-derived negative reference was not created'
    content = neg_ref_files[0].read_text()
    assert '>UnrelatedGene' in content
    for excluded in ('GeneA', 'GeneA_AS', 'GeneA_paralog'):
        assert f'>{excluded}\n' not in content, (
            f'{excluded} was supposed to be excluded from the negative reference'
        )
