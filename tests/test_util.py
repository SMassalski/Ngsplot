import pytest
from mgplot.util import group_by_chromosome, keep_region, filter_regions
from mgplot.core import SignalSegment, Region


def test_group_by_chromosome_regions():
    r1, r2, r3 = Region('chr1', 1, 2, 'Region1', True, 1),\
                 Region('chr1', 1, 2, 'Region2', True, 1),\
                 Region('chr2', 1, 2, 'Region3', True, 1)
    regs = [r1, r2, r3]
    grouped = group_by_chromosome(regs)
    assert r1 in grouped['chr1']
    assert r2 in grouped['chr1']
    assert r3 in grouped['chr2']


def test_group_by_chromosome_signal_segments():
    s1, s2, s3 = SignalSegment('chr1', 1, 2, 1),\
                 SignalSegment('chr1', 1, 2, 1),\
                 SignalSegment('chr2', 1, 2, 1)
    regs = [s1, s2, s3]
    grouped = group_by_chromosome(regs)
    assert s1 in grouped['chr1']
    assert s2 in grouped['chr1']
    assert s3 in grouped['chr2']
    

@pytest.fixture
def sample_region():
    return Region('chr1', 5, 10, 'Region1', True, 10.2)


def test_keep_region_default(sample_region):
    assert keep_region(sample_region) is True


def test_keep_region_omit_chr(sample_region):
    assert keep_region(sample_region, omit_chr=['chr2']) is True
    assert keep_region(sample_region, omit_chr=['chr1']) is False


def test_keep_region_only_chr(sample_region):
    assert keep_region(sample_region, only_chr=['chr2']) is False
    assert keep_region(sample_region, only_chr=['chr1']) is True


def test_keep_region_omit_reg(sample_region):
    assert keep_region(sample_region, omit_reg=['NotRegion1']) is True
    assert keep_region(sample_region, omit_reg=['Region1']) is False


def test_keep_region_min_score(sample_region):
    assert keep_region(sample_region, min_score=10.0) is False
    assert keep_region(sample_region, min_score=10.2) is True
    assert keep_region(sample_region, min_score=10.3) is True


def test_keep_region_max_score(sample_region):
    assert keep_region(sample_region, max_score=10.0) is True
    assert keep_region(sample_region, max_score=10.2) is True
    assert keep_region(sample_region, max_score=10.3) is False


@pytest.fixture
def region_list():
    return [
        Region('chr1', 1, 2, 'Region1', True, 1),
        Region('chr1', 2, 3, 'Region1', True, 1),
        Region('chr4', 1, 2, 'Region2', True, 2),
        Region('chr2', 1, 2, 'Region3', True, 1),
        Region('chr1', 3, 4, 'Region1', True, 2),
        Region('chr2', 1, 2, 'Region4', True, 3),
        Region('chr2', 1, 2, None, True, 0),
        Region('chr3', 4, 9, None, True, 0),
    ]


def test_filter_regions_default(region_list):
    assert region_list == filter_regions(region_list)


def test_filter_regions_only_first(region_list):
    expected = [
        Region('chr1', 1, 2, 'Region1', True, 1),
        Region('chr4', 1, 2, 'Region2', True, 2),
        Region('chr2', 1, 2, 'Region3', True, 1),
        Region('chr2', 1, 2, 'Region4', True, 3),
        Region('chr2', 1, 2, None, True, 0),
        Region('chr3', 4, 9, None, True, 0),
    ]
    assert expected == \
           filter_regions(region_list, only_first=True)


def test_filter_regions_only_first_none(region_list):
    filtered = filter_regions(region_list, only_first=True)
    assert Region('chr2', 1, 2, None, True, 0) in filtered
    assert Region('chr3', 4, 9, None, True, 0) in filtered


def test_filter_regions_n_best_equals_length(region_list):
    assert region_list == filter_regions(region_list, n_best=len(region_list))


def test_filter_regions_n_best_shorter_list(region_list):
    assert region_list == filter_regions(region_list,
                                         n_best=len(region_list)+1)


def test_filter_regions_n_best(region_list):
    expected = [
        Region('chr4', 1, 2, 'Region2', True, 2),
        Region('chr1', 3, 4, 'Region1', True, 2),
        Region('chr2', 1, 2, 'Region4', True, 3),
    ]
    assert expected == filter_regions(region_list, n_best=3)
