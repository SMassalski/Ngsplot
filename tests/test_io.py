import logging
from os.path import abspath, dirname, join
from mgplot.io import BedParser


class TestBedParsers:
    fp = abspath(join(dirname(__file__), 'data', 'hg38_short.testbed'))
    
    def test_parse_read_all(self):
        parser = BedParser()
        results = parser.parse(self.fp)
        assert len(results) == 30
        
    def test_negative_strand_flip(self):
        parser = BedParser()
        results = parser.parse(self.fp)
        assert results[2].start == 29570
        assert results[2].end == 14404
        
    def test_bed_parse_ignore_track_line(self):
        fp = abspath(join(dirname(__file__), 'data',
                          'hg38_short_with_track_line.testbed'))
        parser = BedParser()
        results = parser.parse(fp)
        assert len(results) == 8
        
    def test_bed_parse_missing_field_info(self, caplog):
        caplog.set_level(logging.INFO)
        fp = abspath(join(dirname(__file__), 'data',
                          'hg38_short_with_missing_fields.testbed'))
        parser = BedParser()
        _ = parser.parse(fp)
        assert 'Some optional fields were missing in the file' in caplog.text
        
    def test_parse_region_types(self):
        parser = BedParser()
        results = parser.parse(self.fp)
        region = results[0]
        assert isinstance(region.chromosome, str)
        assert isinstance(region.start, int)
        assert isinstance(region.end, int)
        assert isinstance(region.name, str)
        assert isinstance(region.score, float)
        assert isinstance(region.positive_strand, bool)
