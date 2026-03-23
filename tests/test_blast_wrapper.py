import pytest
import inspect
from hcr_prober import blast_wrapper

def test_best_coverage_uses_args_bitscore():
    """Verify best-coverage strategy uses args.min_bitscore, not hardcoded 75.0."""
    source = inspect.getsource(blast_wrapper.filter_probes_by_blast)
    assert '75.0' not in source, "Hardcoded 75.0 still present in filter_probes_by_blast"

def test_no_bare_exit():
    """Verify no bare exit(1) calls in blast_wrapper."""
    source = inspect.getsource(blast_wrapper)
    # Check for exit(1) that isn't sys.exit(1)
    import re
    bare_exits = re.findall(r'(?<!sys\.)exit\(1\)', source)
    assert len(bare_exits) == 0, f"Found {len(bare_exits)} bare exit(1) calls"

def test_run_negative_screen_signature():
    """Verify run_negative_screen function exists with expected signature."""
    from hcr_prober.blast_wrapper import run_negative_screen
    import inspect
    sig = inspect.signature(run_negative_screen)
    params = list(sig.parameters.keys())
    assert 'probes' in params
    assert 'args' in params
    assert 'temp_dir' in params
    assert 'target_ids' in params
