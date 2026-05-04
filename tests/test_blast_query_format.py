"""Tests for the joined-with-N BLAST query format and strand restriction.

Audit C3 / C5 — keep the joined-arm query (HCR signal requires both arms
to colocalise, so screening for off-targets where both arms hit
adjacently is biologically appropriate), but tighten:

  - Use uppercase 'NN' as the 2-nt gap placeholder (BLAST conventionally
    treats uppercase ambiguity bases). Lowercase is sometimes treated
    as soft-masked depending on -dust settings.
  - Restrict to the plus strand. Probes are antisense to the mRNA, and
    the reference is sense; allowing both strands invites spurious hits
    when the reference contains residual antisense contigs.
"""
import os
from types import SimpleNamespace

from hcr_prober.blast_wrapper import _run_blast


class _StubResult:
    returncode = 0
    stderr = ''


def test_blast_query_uses_uppercase_N_gap_and_plus_strand(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, **kw):
        captured['cmd'] = list(cmd)
        # Find the -query path in the command and read it before the
        # _run_blast 'finally' block deletes it.
        try:
            qi = cmd.index('-query')
            with open(cmd[qi + 1]) as fh:
                captured['query_fasta'] = fh.read()
        except (ValueError, FileNotFoundError):
            captured['query_fasta'] = ''
        # Touch the output path so _run_blast doesn't bail.
        oi = cmd.index('-out')
        open(cmd[oi + 1], 'w').close()
        return _StubResult()

    monkeypatch.setattr('subprocess.run', fake_run)
    probes = [{
        'pair_id': 'test_pair_1',
        'probe_dn_target': 'A' * 25,
        'probe_up_target': 'T' * 25,
    }]
    out = _run_blast(probes, 'fake_db', str(tmp_path), [])
    assert out is not None
    cmd = captured['cmd']
    assert '-strand' in cmd, f'BLAST cmd missing -strand flag: {cmd}'
    assert cmd[cmd.index('-strand') + 1] == 'plus', (
        f'BLAST -strand value should be plus, got: {cmd[cmd.index("-strand") + 1]}'
    )
    fasta = captured['query_fasta']
    assert 'NN' in fasta, f'Joined-arm gap should be uppercase NN; FASTA was:\n{fasta}'
    assert 'nn' not in fasta, (
        f'Lowercase nn should not appear in BLAST query (use uppercase NN). FASTA:\n{fasta}'
    )


def test_negative_screen_uses_real_db_path_not_directory(tmp_path, monkeypatch):
    """run_negative_screen used to pass a directory path (not the actual db
    name prefix returned by create_blast_db) to _run_blast, causing BLAST
    to fail with 'memory map file error' and the screen to silently skip.

    This regression test asserts that the BLAST -db argument used by the
    negative screen is a non-directory path (the actual DB name prefix).
    """
    import os
    from hcr_prober.blast_wrapper import run_negative_screen
    from types import SimpleNamespace

    captured_db_args = []

    def fake_run(cmd, **kw):
        if 'blastn' in cmd[0]:
            try:
                idx = cmd.index('-db')
                captured_db_args.append(cmd[idx + 1])
            except ValueError:
                pass
            oi = cmd.index('-out')
            open(cmd[oi + 1], 'w').close()
        return SimpleNamespace(returncode=0, stderr='', stdout='')

    monkeypatch.setattr('subprocess.run', fake_run)

    # Build a minimal negative reference FASTA and target_ids that don't
    # match any of its sequence ids.
    neg_ref_fasta = tmp_path / 'neg_ref.fasta'
    neg_ref_fasta.write_text('>off_target_1\nACGT' * 30 + '\n')
    args = SimpleNamespace(
        blast_ref=str(neg_ref_fasta),
        blast_negative_ref=None,
        negative_bitscore=None,
        negative_evalue=None,
        min_bitscore=75.0,
        max_evalue=1e-10,
        blast_extra_args=[],
    )
    probes = [{
        'pair_id': 'p1',
        'probe_dn_target': 'A' * 25,
        'probe_up_target': 'T' * 25,
    }]
    run_negative_screen(probes, args, str(tmp_path), {'target_id_not_in_ref'})

    assert captured_db_args, 'no BLAST invocation captured for negative screen'
    db_arg = captured_db_args[-1]
    # The DB arg passed to BLAST must NOT be a bare directory; it must be
    # a path prefix whose basename can be appended to ".nhr"/".nsq" etc.
    assert not os.path.isdir(db_arg), (
        f'Negative-screen BLAST -db points at a directory ({db_arg}); '
        f'should be the DB name prefix returned by create_blast_db.'
    )
