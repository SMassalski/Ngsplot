from collections import namedtuple
from abc import ABC, abstractmethod
import logging
from .util import filter_regions


# DESIGN: The custom parser probably should be in the main script
Region = namedtuple('Region', ['chromosome', 'start', 'end', 'name',
                               'positive_strand', 'score'])


def load_regions(fp, parser='bed', region_part='body', **kwargs):
    r"""Load regions from file

    Parameters
    ----------
    fp : str
        Path of the file containing genomic regions.
    parser : str, one of ['bed', 'score_tsv'] or RegionParserBase
             subclass object
        Parser for the file.
    region_part : str
        Part of the region that will be plotted.
    \*\*kwargs : dict
        Keyword arguments used for filtering. See documentation for
        `filter_regions` and `_keep_region`

    Returns
    -------
    Dict[str,List[Tuple[int,int]]]
        Mapping of chromosomes to starts and ends of regions

    Raises
    ------
    ValueError
        If the provided `parser` str is not recognized or not an object
        subclassing RegionParserBase.
    """
    if type(parser) == str:
        if parser == 'bed':
            parser = BedParser()
            
        elif parser == 'score_tsv':
            parser = ScoreTSVParser()
        else:
            raise ValueError(f'Parser {parser} not recognized.')
    elif not isinstance(parser, RegionParserBase):
        raise ValueError('Invalid parser.')
    
    regions = parser.parse(fp)
    filter_regions(regions, **kwargs)
    # keeps only start and end of regions
    c = {}
    for r in regions:
        c[r.chromosome] = c.get(r.chromosome, []) + [(r.start, r.end)]
    
    # sort by roi
    for key in c:
        if region_part == 'body' or region_part == 'tss':
            c[key] = sorted(c[key])
        else:
            c[key] = sorted(c[key], key=lambda x: x[1])
    return c


class RegionParserBase(ABC):
    """Base class for region parser

    """
    
    @abstractmethod
    def parse(self, fp):
        """Read and parse a file containing genomic regions

        Parameters
        ----------
        fp : str
            Path of the file.

        Returns
        -------
        List[Regions]
            List of regions in the file
        """
        pass


class BedParser(RegionParserBase):
    """Parser for bed files.

    The parser reads optional fields. If not present in the file defaults
    will be used:
    - name: None
    - score: 0
    - strand: +
    """
    
    # TODO: Check if flipping strands make sense
    def parse(self, fp):
        result = []
        log_missing_info = False
        
        with open(fp) as f:
            for line in f:
                line = line.split()
                if line[0].lower() == 'track':  # ignoring track lines
                    continue
                chromosome, start, end = line[:3]
                start = int(start)
                end = int(end)
                name = None
                score = 0
                pos_strand = True
                
                try:
                    name = line[3]
                    score = float(line[4])
                    pos_strand = line[5] == '+'
                except IndexError:
                    log_missing_info = True
                
                # flipping negative strands
                if pos_strand:
                    region = Region(chromosome, start, end, name, pos_strand,
                                    score)
                else:
                    region = Region(chromosome, end, start, name, pos_strand,
                                    score)
                
                result.append(region)
        
        if log_missing_info:
            logging.info('Some optional fields were missing in the file')
        
        return result


class ScoreTSVParser(RegionParserBase):
    """Parser for a custom tsv format"""
    def parse(self, fp):
        result = []
        
        with open(fp, 'r') as f:
            for line in f:
                line = line.split('_')
                name = line[1]
                line = line[-1].split(':')
                chromosome = line[0][:-1]
                line = line[1].split('\t')
                score = int(line[1])
                line = line[0].split('-')
                start = int(line[0])
                end = int(line[1])
                
                result.append(
                    Region(chromosome, start, end, name, True, score))
        
        return result
