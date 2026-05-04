"""Tests for input FASTA validation in file_io.read_fasta.

Pre-fix: an empty FASTA returned an empty dict silently and the caller
in main.py would crash later with a confusing error. Malformed FASTA
fell through to the same silent behaviour.

After Task 2.2: empty input raises ValueError; the caller can produce
a clear error early.
"""
import pytest


def test_empty_fasta_raises(tmp_path):
    fa = tmp_path / 'empty.fasta'
    fa.write_text('')
    from hcr_prober.file_io import read_fasta
    with pytest.raises(ValueError, match='no sequences'):
        read_fasta(str(fa))


def test_fasta_with_only_header_no_sequence_raises(tmp_path):
    """A FASTA with a header but no sequence body is structurally degenerate."""
    fa = tmp_path / 'header_only.fasta'
    fa.write_text('>my_seq\n')
    from hcr_prober.file_io import read_fasta
    with pytest.raises(ValueError):
        read_fasta(str(fa))


def test_valid_fasta_round_trips(tmp_path):
    fa = tmp_path / 'good.fasta'
    fa.write_text('>seq1\n' + 'ACGT' * 30 + '\n>seq2\n' + 'TGCA' * 30 + '\n')
    from hcr_prober.file_io import read_fasta
    seqs = read_fasta(str(fa))
    assert set(seqs.keys()) == {'seq1', 'seq2'}
    assert len(seqs['seq1']) == 4 * 30
