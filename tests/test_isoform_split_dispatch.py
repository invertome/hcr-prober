"""TDD test for v1.10.0 audit finding C1.

Bug: ``main.py`` references ``args._isoform_ids`` in ``create_probe_blueprint``
to extend the negative-BLAST target_ids set, but the isoform-split dispatch
never actually sets that attribute. Result: probes for an isoform's ``common``
region get rejected by the auto-derived negative screen because the screen
treats sister isoforms as off-targets.

This test runs the isoform-split command path with everything heavy mocked
out, capturing the args object passed to ``create_probe_blueprint``. After
the fix, every captured args (both the common job and each unique job)
must have ``_isoform_ids`` populated.
"""
import sys
import pytest
from types import SimpleNamespace


COMMON_REGION = (
    'ATGCAGCATGGAGTCGAATCAGCTCAGCATCAGCTAGCATCGATCGATCGAATCGATCGAATG'
    'CGCATGCATGCAGCATCGATCGAGCAGCATCGAGCATCGAGCATCGTAGCATCGATCGAGCAT'
    'CGAATCGATCGATCGAATCGATCAGCATCGATCGATCAGCATCGATCAGCATGCATCAGCAGC'
    'ATCGATCGATCAGCATGCAGCATGCAGCATGCATGCATCGATCAGCATCGAATCGAATCGCGT'
)


@pytest.fixture
def isoform_fasta(tmp_path):
    """Three isoforms of `gene1` sharing a >200 bp common 5' region."""
    fa = tmp_path / 'isoforms.fasta'
    fa.write_text(
        f'>gene1_iso1\n{COMMON_REGION}AAAAATTTTTGGGGGCCCCC' * 1 + 'AAAAATTTTT' * 5 + '\n'
        f'>gene1_iso2\n{COMMON_REGION}TTTTTGGGGGCCCCCAAAAA' * 1 + 'TTTTTGGGGG' * 5 + '\n'
        f'>gene1_iso3\n{COMMON_REGION}GGGGGCCCCCAAAAATTTTT' * 1 + 'GGGGGCCCCC' * 5 + '\n'
    )
    return fa


def test_isoform_split_sets_isoform_ids_on_dispatched_args(isoform_fasta, tmp_path, monkeypatch):
    captured = []

    def mock_create_blueprint(gene_name, seq, td, args):
        captured.append({
            'gene_name': gene_name,
            'isoform_ids': getattr(args, '_isoform_ids', None),
        })
        return None, None, {}

    monkeypatch.setattr('hcr_prober.main.create_probe_blueprint', mock_create_blueprint)
    monkeypatch.setattr('hcr_prober.main.check_dependencies', lambda: None)
    monkeypatch.setattr('hcr_prober.blast_wrapper.create_blast_db', lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.prober.finalize_probes', lambda *a, **kw: [])
    monkeypatch.setattr('hcr_prober.prober.subsample_probes', lambda probes, n: probes)
    monkeypatch.setattr('hcr_prober.file_io.write_outputs', lambda *a, **kw: None)
    monkeypatch.setattr(sys, 'argv', [
        'hcr-prober', 'isoform-split',
        '-i', str(isoform_fasta),
        '-o', str(tmp_path / 'out'),
        '--gene-prefix', 'gene1',
        '--amplifier', 'B1',
    ])

    from hcr_prober.main import main
    main()

    assert captured, 'create_probe_blueprint was never invoked'

    common_jobs = [c for c in captured if c['gene_name'].endswith('_common')]
    unique_jobs = [c for c in captured if not c['gene_name'].endswith('_common')]

    assert common_jobs, 'No common-region job was dispatched'
    assert unique_jobs, 'No unique-region jobs were dispatched'

    expected_iso_set = {'gene1_iso1', 'gene1_iso2', 'gene1_iso3'}

    for job in common_jobs:
        assert job['isoform_ids'] == expected_iso_set, (
            f"Common job '{job['gene_name']}' must include all sister isoforms in "
            f"_isoform_ids so the negative screen does not reject probes that "
            f"legitimately bind them.\n  expected: {expected_iso_set}\n  got: {job['isoform_ids']}"
        )

    for job in unique_jobs:
        expected_for_unique = expected_iso_set - {job['gene_name']}
        assert job['isoform_ids'] == expected_for_unique, (
            f"Unique job '{job['gene_name']}' must exclude its OWN id from _isoform_ids "
            f"(otherwise the negative screen would treat self-hits as off-target).\n"
            f"  expected: {expected_for_unique}\n  got: {job['isoform_ids']}"
        )
