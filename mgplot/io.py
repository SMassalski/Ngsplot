from collections import namedtuple
from abc import ABC, abstractmethod
import logging


# DESIGN: Filtering does not belong in this module
# DESIGN: The custom parser probably should be in the main script
Region = namedtuple('Region', ['chromosome', 'start', 'end', 'name',
                               'positive_strand', 'score'])


def load_regions(fp, parser='bed', region_part='body', **kwargs):
    r"""Load regions from file

    Parameters
    ----------
    fp : str
        Path of the bed file.
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
    Dict[str:List[Tuple[int,int]]]
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


def _keep_region(region, omit_chr=(), omit_reg=(), only_chr=None,
                 max_score=None, min_score=None, **_):
    r"""Check if a region satisfies filtering conditions.

    Parameters
    ----------
    region : Region
        Region to be checked.
    omit_chr : Iterable[str]
        An iterable containing names of chromosomes to be omitted.
    omit_reg : Iterable[str]
        An iterable containing names of regions to be omitted. This
        argument will be ignored if the bed file does not contain region
        names.
    only_chr : Iterable[str]
        An iterable of chromosome names that will be exclusively
         considered.
    max_score : float
        Upper limit (inclusive) of scores. Only regions with scores
        lower or equal will be returned.
    min_score : float
        Lower limit (inclusive) of scores. Only regions with scores
        higher or equal will be returned.

    Returns
    -------
    bool
        Whether the region passes the check.
    """
    # selecting chromosomes
    if only_chr is not None and region.chromosome not in only_chr:
        return False
    
    # omitting chromosomes
    if region.chromosome in omit_chr:
        return False
    
    # omitting regions
    if region.name in omit_reg:
        return False
    
    # score range filter
    if max_score is not None and region.score <= max_score:
        return False
    if min_score is not None and region.score >= min_score:
        return False
    
    return True


def filter_regions(regions, only_first=False,
                   n_best=None, **kwargs):
    r"""Filter list of genomic regions.

    Parameters
    ----------
    regions : List[Region]
        Regions to be filtered.
    
    only_first : bool
        Whether to keep only the first region when the bed file
        contains multiple regions with the same name. This argument
        will be ignored if the bed file does not contain region
        names.
    n_best : int
        Number of regions with the highest score to be kept.
    \*\*kwargs : dict
        Additional keyword arguments including elementwise filtering
        conditions. Valid elementwise conditions include:
            - omit_chr : Iterable[str]
            - omit_reg : Iterable[str]
            - only_chr : Iterable[str]
            - max_score : float
            - min_score : float
        See documentation for _keep_region

    Returns
    -------
    List[Region]
        Filtered list of regions.
    """
    # elementwise filter
    regions = list(filter(lambda x: _keep_region(x, **kwargs), regions))
    result = []
    
    # only_first filter
    seen_names = set()
    if only_first:
        for r in regions:
            if r.name in seen_names:
                continue
            else:
                if r.name is not None:
                    seen_names.add(r.name)
                result.append(r)
    
    # n_best filter
    if n_best is not None:
        if n_best < len(result):
            # keeping the previous order
            n_best_regions = sorted(result,
                                    key=lambda t: t.score,
                                    reverse=True)
            n_best_regions = set(n_best_regions[:n_best])
            result = list(filter(lambda x: x in n_best_regions, result))
    return result
