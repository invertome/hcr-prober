import pytest
from hcr_prober.utils.thermo_utils import calculate_gc_content, calculate_tm

def test_gc_content_basic():
    assert calculate_gc_content("GCGCGC") == pytest.approx(100.0)
    assert calculate_gc_content("ATATAT") == pytest.approx(0.0)
    assert calculate_gc_content("GCATAT") == pytest.approx(33.33, rel=0.01)

def test_gc_content_empty():
    assert calculate_gc_content("") == 0.0

def test_tm_returns_float():
    tm = calculate_tm("ATGCGATCGATCGATCGATCGATCG")
    assert isinstance(tm, float)
    assert 40.0 < tm < 80.0

def test_tm_empty():
    assert calculate_tm("") == 0.0

def test_tm_custom_concentrations():
    tm_default = calculate_tm("ATGCGATCGATCGATCGATCGATCG")
    tm_high_na = calculate_tm("ATGCGATCGATCGATCGATCGATCG", Na=200)
    # Higher Na+ should increase Tm
    assert tm_high_na > tm_default

def test_tm_caching():
    """Calling with same args should return cached result."""
    seq = "ATGCGATCGATCGATCGATCGATCG"
    r1 = calculate_tm(seq)
    r2 = calculate_tm(seq)
    assert r1 == r2
    # Verify cache is working by checking cache_info
    info = calculate_tm.cache_info()
    assert info.hits >= 1
