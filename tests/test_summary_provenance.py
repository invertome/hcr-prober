"""Tests for the provenance section in the run summary file.

Phase 4.1. The summary needs to embed enough context that any single
output file is reproducible. That means: the exact command line invoked,
the version of hcr-prober and its critical dependencies, the host the
run was made on, when it ran, and what RNG seed was in effect.
"""
import os
import re
import shutil
import subprocess

import pytest


pytestmark = pytest.mark.skipif(
    not shutil.which('hcr-prober'),
    reason='hcr-prober console script not on PATH',
)


def _run(tmp_path):
    target = os.path.join(os.path.dirname(__file__), 'fixtures', 'mock_target.fasta')
    out_dir = tmp_path / 'out'
    cmd = [
        'hcr-prober', 'design',
        '-i', target,
        '-o', str(out_dir),
        '--amplifier', 'B1',
        '--max-probes', '3',
        '--seed', '7',
        '--min-gc', '20', '--max-gc', '80',
        '--min-tm', '0', '--max-tm', '100',
        '--skip-5prime', '0',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result, out_dir


def _summary_text(out_dir):
    for root, _, files in os.walk(out_dir):
        for f in files:
            if f.endswith('_summary.txt'):
                with open(os.path.join(root, f)) as fh:
                    return fh.read()
    pytest.fail(f'No summary file under {out_dir}')


def test_summary_contains_command_line(tmp_path):
    result, out_dir = _run(tmp_path)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'Command:' in text or 'argv' in text.lower()
    assert '--amplifier' in text
    assert '--seed' in text and '7' in text


def test_summary_contains_versions_and_seed(tmp_path):
    result, out_dir = _run(tmp_path)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert 'hcr-prober' in text
    assert re.search(r'\b\d+\.\d+\.\d+\b', text), 'no semver in summary'
    assert re.search(r'BLAST\+?[^\n]*\d+\.\d+', text), 'BLAST+ version missing in summary'
    assert re.search(r'primer3[^\n]*\d', text), 'primer3 version missing in summary'
    assert re.search(r'[Ss]eed[^\n]*7', text), '--seed value missing in summary'


def test_summary_contains_timestamp_and_host(tmp_path):
    result, out_dir = _run(tmp_path)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    text = _summary_text(out_dir)
    assert re.search(r'\d{4}-\d{2}-\d{2}', text), 'date stamp missing in summary'
    assert 'host' in text.lower() or 'Host' in text, 'hostname missing in summary'
