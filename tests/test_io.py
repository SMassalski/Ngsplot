import logging
import os

import pytest

from mgplot.io import BedParser, read_bedgraph
from mgplot.core import SignalSegment


data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


class TestBedParsers:
    
    main_fp = os.path.join(data_dir, 'hg38_short.testbed')
    track_line_fp = os.path.join(data_dir,
                                 'hg38_short_with_track_line.testbed')
    missing_field_fp = os.path.join(data_dir,
                                    'hg38_short_with_missing_fields.testbed')
    
    @pytest.fixture(scope='class', autouse=True)
    def parser(self):
        return BedParser()
    
    def test_parse_read_all(self, parser):
        results = parser.parse(self.main_fp)
        assert len(results) == 30
    
    def test_bed_parse_ignore_track_line(self, parser):
        results = parser.parse(self.track_line_fp)
        assert len(results) == 8
        
    def test_bed_parse_missing_field_info(self, caplog, parser):
        caplog.set_level(logging.INFO)
        _ = parser.parse(self.missing_field_fp)
        assert 'Some optional fields were missing in the file' in caplog.text
        
    def test_parse_region_types(self, parser):
        results = parser.parse(self.main_fp)
        region = results[0]
        assert isinstance(region.chromosome, str)
        assert isinstance(region.start, int)
        assert isinstance(region.end, int)
        assert isinstance(region.name, str)
        assert isinstance(region.score, float)
        assert isinstance(region.positive_strand, bool)


@pytest.fixture
def bedgraph_segments():
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


def test_read_bedgraph(bedgraph_segments):
    assert read_bedgraph(os.path.join(data_dir, 'test.bedGraph')) == \
           bedgraph_segments


def test_read_bedgraph_with_track_lines(bedgraph_segments):
    read = read_bedgraph(os.path.join(data_dir,
                                      'test_with_track_lines.bedGraph'))
    assert read == bedgraph_segments
