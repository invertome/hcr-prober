"""Tests for the BLAST DB cache key.

Audit C7. Pre-fix: create_blast_db used os.path.splitext(os.path.basename(...))
as the DB name. Two different reference FASTAs with the same basename in
two different directories would collide in a shared --db-path, silently
overwriting one DB with the other.

After the fix: a content-derived component (e.g., short sha256 of the
absolute ref-fasta path) is mixed into the DB name so the cache key is
unique per source path.
"""
import os
import shutil
import pytest


def test_two_refs_with_same_basename_get_distinct_db_names(tmp_path, monkeypatch):
    """Two FASTA files in different directories sharing the same basename
    must yield distinct DB names so they do not collide in --db-path."""
    # Avoid actually invoking makeblastdb during the test.
    monkeypatch.setattr('subprocess.run', lambda *a, **kw: type('R', (), {'returncode': 0, 'stderr': '', 'stdout': ''})())

    from hcr_prober.blast_wrapper import create_blast_db

    dir_a = tmp_path / 'a'
    dir_b = tmp_path / 'b'
    dir_a.mkdir()
    dir_b.mkdir()
    (dir_a / 'transcriptome.fasta').write_text('>seq_a\nACGT\n' * 5)
    (dir_b / 'transcriptome.fasta').write_text('>seq_b\nTGCA\n' * 5)

    db_path = tmp_path / 'shared_blast_dbs'

    name_a = create_blast_db(str(dir_a / 'transcriptome.fasta'), str(db_path))
    name_b = create_blast_db(str(dir_b / 'transcriptome.fasta'), str(db_path))

    assert name_a != name_b, (
        f'Two refs with same basename collided to identical DB name: {name_a!r}. '
        f'Cache key must include path-derived component.'
    )
