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
