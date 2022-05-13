import logging
import os

import pytest

from mgplot.io import BedParser


class TestBedParsers:
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
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
    
    @pytest.mark.skip(reason='Handling of negative strands undecided')
    def test_negative_strand_flip(self, parser):
        results = parser.parse(self.main_fp)
        assert results[2].start == 29570
        assert results[2].end == 14404

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
