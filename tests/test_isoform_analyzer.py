import pytest
from hcr_prober.isoform_analyzer import find_common_regions

def test_find_common_regions_identical():
    """Identical sequences should have one common region spanning most of the length."""
    seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    group = {'iso1': seq, 'iso2': seq}
    ref_id, intervals = find_common_regions(group)
    assert len(intervals) > 0
    total_covered = sum(end - start for start, end in intervals)
    assert total_covered > len(seq) * 0.5

def test_find_common_regions_no_similarity():
    """Completely different sequences should have no common regions."""
    group = {
        'iso1': 'A' * 200,
        'iso2': 'G' * 200 + 'C' * 200,
    }
    ref_id, intervals = find_common_regions(group)
    total_covered = sum(end - start for start, end in intervals)
    assert total_covered < 52

def test_no_deprecation_import():
    """Verify Bio.pairwise2 is not imported."""
    import inspect
    import hcr_prober.isoform_analyzer as ia
    source = inspect.getsource(ia)
    assert 'pairwise2' not in source, "Bio.pairwise2 still referenced"

def test_single_isoform():
    """Single isoform should return full sequence as common."""
    group = {'iso1': 'ATGCGATCGATCGATCGATCGATCG' * 5}
    ref_id, intervals = find_common_regions(group)
    assert ref_id == 'iso1'
    assert len(intervals) == 1
    assert intervals[0] == (0, len(group['iso1']))

def test_min_interval_len_parameterized():
    """Verify min_interval_len parameter is accepted."""
    seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    group = {'iso1': seq, 'iso2': seq}
    # Should work with custom min_interval_len
    ref_id, intervals = find_common_regions(group, min_interval_len=10)
    assert len(intervals) > 0
