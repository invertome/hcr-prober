"""End-to-end integration tests exercising the BLAST screen.

Phase 3.1. The unit tests in this repo mock subprocess.run and never
actually invoke blastn/makeblastdb. These tests run the real pipeline
against the mock_transcriptome.fasta / mock_target.fasta fixtures and
exercise each --positive-selection-strategy plus the negative-screen
auto-derivation path.
"""
import os
import shutil
import subprocess

import pytest


FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')
TARGET_FASTA = os.path.join(FIXTURE_DIR, 'mock_target.fasta')
TRANSCRIPTOME_FASTA = os.path.join(FIXTURE_DIR, 'mock_transcriptome.fasta')

pytestmark = pytest.mark.skipif(
    not (shutil.which('blastn') and shutil.which('makeblastdb')),
    reason='NCBI BLAST+ not on PATH; integration tests require it',
)


def _run(tmp_path, *extra):
    out_dir = tmp_path / 'out'
    db_dir = tmp_path / 'blast_db'
    cmd = [
        'hcr-prober', 'design',
        '-i', TARGET_FASTA,
        '-o', str(out_dir),
        '--db-path', str(db_dir),
        '--blast-ref', TRANSCRIPTOME_FASTA,
        '--amplifier', 'B1',
        '--max-probes', '5',
        '--seed', '42',
        '--min-gc', '20', '--max-gc', '80',
        '--min-tm', '0', '--max-tm', '100',
        '--skip-5prime', '0',
    ] + list(extra)
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, out_dir


def _summary_text(out_dir):
    for root, _, files in os.walk(out_dir):
        for f in files:
            if f.endswith('_summary.txt'):
                with open(os.path.join(root, f)) as fh:
                    return fh.read()
    pytest.fail(f'No summary file under {out_dir}')


def test_blast_any_strong_hit_strategy(tmp_path):
    result, out_dir = _run(tmp_path, '--positive-selection-strategy', 'any-strong-hit')
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'Specific Candidates (Post-BLAST)' in text
    assert 'any-strong-hit' in text


def test_blast_best_coverage_strategy(tmp_path):
    result, out_dir = _run(tmp_path, '--positive-selection-strategy', 'best-coverage')
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'best-coverage' in text


def test_blast_specific_id_strategy(tmp_path):
    """specific-id requires --target-transcript-id; test that it runs and reports the strategy."""
    result, out_dir = _run(
        tmp_path,
        '--positive-selection-strategy', 'specific-id',
        '--target-transcript-id', 'mock_gene_alpha',
    )
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'specific-id' in text


def test_blast_specific_id_without_target_id_exits_with_error(tmp_path):
    """Missing --target-transcript-id should fail with a clear FATAL error."""
    result, _ = _run(tmp_path, '--positive-selection-strategy', 'specific-id')
    assert result.returncode != 0
    assert 'specific-id' in result.stderr.lower()
    assert 'target-transcript-id' in result.stderr.lower()


def test_negative_screen_runs_in_audit_funnel(tmp_path):
    """When --blast-ref is set, the audit funnel should include the negative-BLAST step."""
    result, out_dir = _run(tmp_path)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'After Negative BLAST Screen' in text


def test_blast_ref_with_target_id_excluded_from_negative(tmp_path):
    """When the target sequence is in the BLAST reference, the auto-derived
    negative reference must exclude it. Otherwise every probe gets rejected
    (the original audit-driven scenario)."""
    # mock_target.fasta is mock_gene_alpha; mock_transcriptome.fasta also
    # contains mock_gene_alpha. Default behaviour (target_ids = {gene_name})
    # should already cover this.
    result, out_dir = _run(tmp_path)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    # We expect at least the BLAST step values to be non-zero
    import re
    m = re.search(r'(\d+)\s+After Negative BLAST Screen', text)
    assert m, f'Could not find negative-screen line in summary:\n{text}'
    # Need at least one probe to survive (mock_gene_beta and mock_gene_gamma
    # are unrelated random sequences with no significant homology to alpha).
    assert int(m.group(1)) > 0, (
        f'Negative screen rejected all probes — auto-exclude of target may be broken.\n'
        f'Summary excerpt:\n{text}'
    )
