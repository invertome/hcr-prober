"""End-to-end test for isoform-split with real BLAST.

Phase 3.2. Augments the unit test in test_isoform_split_dispatch (which
mocks create_probe_blueprint to confirm _isoform_ids is set) with an
actual subprocess invocation that runs the full pipeline including BLAST.

The fixture: 3 isoforms of one synthetic gene sharing a long 5' common
region, plus 1-2 unique 3' tails per isoform. The BLAST reference is the
isoform set itself plus an unrelated transcript. With the C1 fix in
place, common-region probes that hit sister isoforms must NOT be
rejected by the auto-derived negative screen — the C1 bug would have
removed the entire common-probe set.
"""
import os
import shutil
import subprocess

import pytest


pytestmark = pytest.mark.skipif(
    not (shutil.which('blastn') and shutil.which('makeblastdb')),
    reason='NCBI BLAST+ not on PATH; integration tests require it',
)


COMMON = (
    'ATGCAGCATGGAGTCGAATCAGCTCAGCATCAGCTAGCATCGATCGATCGAATCGATCGAATG'
    'CGCATGCATGCAGCATCGATCGAGCAGCATCGAGCATCGAGCATCGTAGCATCGATCGAGCAT'
    'CGAATCGATCGATCGAATCGATCAGCATCGATCGATCAGCATCGATCAGCATGCATCAGCAGC'
    'ATCGATCGATCAGCATGCAGCATGCAGCATGCATGCATCGATCAGCATCGAATCGAATCGCGT'
    'TCGATCGATCGATCAGCATGCATCAGCATCGATCAGCATCGATCAGCATGCATGCATCGATCG'
)  # ~310 nt common region

UNIQUE_TAILS = {
    'gene1_iso1': 'AAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCC' + 'A' * 50,
    'gene1_iso2': 'CCCCGGGGAAAATTTTCCCCGGGGAAAATTTT' + 'G' * 50,
    'gene1_iso3': 'TTTTCCCCGGGGAAAATTTTCCCCGGGGAAAA' + 'C' * 50,
}


@pytest.fixture
def isoform_workspace(tmp_path):
    """Build an input FASTA + matching BLAST transcriptome."""
    input_fa = tmp_path / 'input.fasta'
    transcriptome = tmp_path / 'transcriptome.fasta'

    with open(input_fa, 'w') as fi, open(transcriptome, 'w') as ft:
        for iso_id, tail in UNIQUE_TAILS.items():
            seq = COMMON + tail
            for fh in (fi, ft):
                fh.write(f'>{iso_id}\n{seq}\n')
        ft.write('>unrelated_gene\n' + 'GCATGCATGCAT' * 80 + '\n')

    return {
        'input': str(input_fa),
        'transcriptome': str(transcriptome),
        'out': str(tmp_path / 'out'),
        'db_path': str(tmp_path / 'blast_db'),
    }


def test_isoform_split_runs_end_to_end(isoform_workspace):
    cmd = [
        'hcr-prober', 'isoform-split',
        '-i', isoform_workspace['input'],
        '-o', isoform_workspace['out'],
        '--db-path', isoform_workspace['db_path'],
        '--blast-ref', isoform_workspace['transcriptome'],
        '--gene-prefix', 'gene1',
        '--amplifier', 'B1',
        '--max-probes', '3',
        '--seed', '42',
        '--min-gc', '20', '--max-gc', '80',
        '--skip-5prime', '0',
        '--max-hairpin-dg', '-50', '--max-homodimer-dg', '-50', '--max-heterodimer-dg', '-50',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'
    out = isoform_workspace['out']
    assert os.path.isdir(os.path.join(out, 'gene1', 'common_probes'))
    assert os.path.isdir(os.path.join(out, 'gene1', 'isoform_specific_probes'))


def test_isoform_split_common_probes_survive_negative_screen(isoform_workspace):
    """C1 regression: probes for the common region hit all sister isoforms by
    design. The negative screen must not reject them. _isoform_ids must
    be propagated so all sisters are excluded from the auto-derived
    negative reference.
    """
    cmd = [
        'hcr-prober', 'isoform-split',
        '-i', isoform_workspace['input'],
        '-o', isoform_workspace['out'],
        '--db-path', isoform_workspace['db_path'],
        '--blast-ref', isoform_workspace['transcriptome'],
        '--gene-prefix', 'gene1',
        '--amplifier', 'B1',
        '--max-probes', '3',
        '--seed', '42',
        '--min-gc', '20', '--max-gc', '80',
        '--skip-5prime', '0',
        '--max-hairpin-dg', '-50', '--max-homodimer-dg', '-50', '--max-heterodimer-dg', '-50',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f'stderr:\n{result.stderr}'

    summary_path = os.path.join(
        isoform_workspace['out'], 'gene1', 'common_probes',
        'gene1_common', 'B1', 'gene1_common_B1_summary.txt',
    )
    assert os.path.exists(summary_path), f'common-probe summary missing: {summary_path}'
    with open(summary_path) as fh:
        text = fh.read()

    import re
    m_post = re.search(r'(\d+)\s+After Negative BLAST Screen', text)
    m_pre = re.search(r'(\d+)\s+Specific Candidates \(Post-BLAST\)', text)
    assert m_post and m_pre, f'Audit funnel missing expected lines:\n{text}'
    pre = int(m_pre.group(1))
    post = int(m_post.group(1))
    assert post > 0, (
        f'Negative screen rejected ALL common-region probes ({pre} -> 0). '
        f'C1 fix has regressed: _isoform_ids must be propagated so sister '
        f'isoforms do not look like off-targets.\n\nSummary excerpt:\n{text[-1500:]}'
    )
