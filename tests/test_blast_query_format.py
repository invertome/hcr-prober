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
    assert cmd[cmd.index('-strand') + 1] == 'minus', (
        f'BLAST -strand value should be minus (probe is antisense; on-target '
        f'match is the revcomp of the probe), got: {cmd[cmd.index("-strand") + 1]}'
    )
    fasta = captured['query_fasta']
    assert 'NN' in fasta, f'Joined-arm gap should be uppercase NN; FASTA was:\n{fasta}'
    assert 'nn' not in fasta, (
        f'Lowercase nn should not appear in BLAST query (use uppercase NN). FASTA:\n{fasta}'
    )


def test_blast_query_includes_num_threads_when_given(tmp_path, monkeypatch):
    """Phase 5.1: --threads N must propagate to blastn -num_threads."""
    captured_cmd = []

    def fake_run(cmd, **kw):
        captured_cmd.append(list(cmd))
        oi = cmd.index('-out')
        open(cmd[oi + 1], 'w').close()
        return _StubResult()

    monkeypatch.setattr('subprocess.run', fake_run)
    probes = [{'pair_id': 'x', 'probe_dn_target': 'A' * 25, 'probe_up_target': 'T' * 25}]
    _run_blast(probes, 'fake_db', str(tmp_path), ['-num_threads', '4'])
    cmd = captured_cmd[0]
    assert '-num_threads' in cmd
    assert cmd[cmd.index('-num_threads') + 1] == '4'


def test_threads_flag_propagates_to_blast_extra_args(monkeypatch):
    """End-to-end: parse --threads 8 and verify args.blast_extra_args picks
    up '-num_threads 8'. Mocks the rest of main() so only the parsing path
    runs."""
    import sys
    from types import SimpleNamespace

    captured_extra = []

    def fake_blueprint(gene_name, seq, td, args):
        captured_extra.append(list(getattr(args, 'blast_extra_args', [])))
        return None, None, {}

    monkeypatch.setattr('hcr_prober.main.create_probe_blueprint', fake_blueprint)
    monkeypatch.setattr('hcr_prober.main.check_dependencies', lambda: None)
    monkeypatch.setattr('hcr_prober.blast_wrapper.create_blast_db',
                        lambda *a, **kw: 'fake_db')
    monkeypatch.setattr('hcr_prober.prober.finalize_probes', lambda *a, **kw: [])
    monkeypatch.setattr('hcr_prober.prober.subsample_probes', lambda probes, n: probes)
    monkeypatch.setattr('hcr_prober.file_io.write_outputs', lambda *a, **kw: None)
    monkeypatch.setattr('hcr_prober.file_io.read_fasta',
                        lambda *a, **kw: {'gene1': 'ACGT' * 100})

    monkeypatch.setattr(sys, 'argv', [
        'hcr-prober', 'design',
        '-i', '/dev/null', '-o', '/tmp/dummy_out',
        '--amplifier', 'B1',
        '--threads', '8',
    ])
    from hcr_prober.main import main
    main()

    assert captured_extra, 'create_probe_blueprint never received args'
    extra = captured_extra[0]
    assert '-num_threads' in extra and '8' in extra, (
        f"expected -num_threads 8 in args.blast_extra_args, got {extra}"
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
