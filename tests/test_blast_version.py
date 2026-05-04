"""Tests for BLAST+ version detection in main.check_dependencies.

Phase 2.5. The dependency check should detect not only that blastn is
present but also which version it is, log the version, and warn if it
falls outside the tested range. This makes incompatibilities surface
early rather than mid-pipeline as cryptic 'Database memory map' or
'-task' errors.
"""
from types import SimpleNamespace
from loguru import logger


def test_check_dependencies_logs_blastn_version(monkeypatch):
    captured_messages = []

    handler_id = logger.add(lambda m: captured_messages.append(m), level='INFO')
    try:
        monkeypatch.setattr('shutil.which', lambda cmd: f'/fake/{cmd}')

        def fake_run(cmd, **kw):
            if cmd[:2] == ['blastn', '-version']:
                return SimpleNamespace(returncode=0, stdout='blastn: 2.16.0+\n Package: blast 2.16.0\n', stderr='')
            return SimpleNamespace(returncode=0, stdout='', stderr='')

        monkeypatch.setattr('subprocess.run', fake_run)

        from hcr_prober.main import check_dependencies
        check_dependencies()
    finally:
        logger.remove(handler_id)

    text = ''.join(captured_messages)
    assert '2.16.0' in text, (
        f'Expected blastn version in log; got:\n{text}'
    )
