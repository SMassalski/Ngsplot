import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

from mgplot.core import normalize, Region, SignalSegment, query


def test_normalize_same_length(random_array):
    assert_allclose(normalize(random_array, len(random_array)),
                    random_array)


@pytest.fixture
def random_array():
    rng = np.random.default_rng()
    return rng.random(50)


@pytest.fixture
def signal():
    return [
        SignalSegment('chr1', 20, 277, 5.16),
        SignalSegment('chr1', 277, 495, 7.92),
        SignalSegment('chr1', 495, 500, -4.22),
        SignalSegment('chr5', 20, 291, -9.35),
        SignalSegment('chr5', 291, 480, 5.74),
        SignalSegment('chr5', 480, 500, 5.99),
        SignalSegment('chrM', 20, 79, 6.04),
        SignalSegment('chrM', 79, 226, -8.45),
        SignalSegment('chrM', 226, 319, 5.0),
        SignalSegment('chrM', 319, 500, -9.99)
    ]


@pytest.fixture
def regions():
    return [
        Region('chr1', 0, 20, 'dummy_region', True, 0),
        Region('chrM', 80, 319, 'dummy_region1', True, 0),
        Region('chr1', 15, 40, 'dummy_region2', False, 0)
    ]


def test_query_defaults(signal):
    region = Region('chr1', 0, 27, 'dummy_region', True, 0)
    result = query([region], signal)[0]
    expected = np.zeros(2001, dtype=result.dtype)
    expected[1020:1277] = 5.16
    expected[1277:1495] = 7.92
    expected[1495:1500] = -4.22
    assert_equal(result, expected)


def test_query_default_body(signal, regions):
    assert query(regions, signal, roi='body', flank=1000).shape == (3, 4000)
    assert query(regions, signal, roi='body', flank=10).shape == (3, 1020)
    assert query(regions, signal, roi='body', flank=0).shape == (3, 1000)


def test_query_start(signal, regions):
    result = query(regions, signal, roi='start', flank=10)
    expected = np.zeros((3, 21), dtype=result.dtype)
    expected[1, 15:] = 5.16
    expected[2, :9] = 6.04
    expected[2, 9:] = -8.45
    assert_equal(result, expected)


def test_query_end(signal, regions):
    result = query(regions, signal, roi='end', flank=10)
    expected = np.zeros((3, 21), dtype=result.dtype)
    expected[0, 10:] = 5.16
    expected[1, :] = 5.16
    expected[2, :10] = 5
    expected[2, 10:] = -9.99
    assert_equal(result, expected)


def test_query_lengths(signal, regions):
    assert query(regions, signal, roi='start', flank=10).shape == (3, 21)
    assert query(regions, signal, roi='end', flank=10).shape == (3, 21)
    assert query(regions, signal, roi='body', flank=10, body=10).shape == \
           (3, 30)


def test_query_flip_strands(signal, regions):
    result = query(regions[2:], signal, flank=10, flip_negative=True)[0]
    expected = np.zeros(21, dtype=result.dtype)
    expected[:-15] = 5.16
    assert_equal(result, expected)


def test_query_incorrect_roi(signal, regions):
    with pytest.raises(ValueError):
        query(regions, signal, roi='not_real_roi')


def test_query_missing_chromosome(signal, regions, caplog):
    regions.append(Region('fake_chr', 15, 40, 'dummy_region2', False, 0))
    query(regions, signal)
    assert 'Signal info missing for chromosome fake_chr' in caplog.text
