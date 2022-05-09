from collections import namedtuple
import logging

# TODO: Implement parser class
# TODO: Separate filtering from parsing

Region = namedtuple('Region', ['chromosome', 'start', 'end', 'name',
                               'positive_strand', 'score'])


def load_regions(fp, file_format='bed', region_part='body', omit_chr=(),
                 omit_reg=(), only_chr=None, only_first=False, n_best=None,
                 max_score=None, min_score=None):
    """Load regions from file

    Parameters
    ----------
    fp : str
        Path of the bed file.
    file_format : str, one of ['bed', 'score_tsv']
        Format of the input file.
    region_part : str
        Part of the region that will be plotted.
    omit_chr : Iterable[str]
        An iterable containing names of chromosomes to be omitted.
    omit_reg : Iterable[str]
        An iterable containing names of regions to be omitted. This
        argument will be ignored if the bed file does not contain region
        names.
    only_first : bool
        Whether to keep only the first region when the bed file contains
        multiple regions with the same name. This argument will be
        ignored if the bed file does not contain region names.
    only_chr : Iterable[str]
        An iterable of chromosome names that will be exclusively
         considered.
    n_best : int
        Number of regions with the highest score to be returned.
    max_score : float
        Upper limit (inclusive) of scores. Only regions with scores
        lower or equal will be returned.
    min_score : float
        Lower limit (inclusive) of scores. Only regions with scores
        higher or equal will be returned.

    Returns
    -------
    Dict[str:List[Tuple[int,int]]]
        Mapping of chromosomes to starts and ends of regions

    Raises
    ------
    ValueError
        If the provided `file_format` is not supported.
    """
    if file_format == 'bed':
        regions = read_bed(fp, omit_chr=omit_chr, omit_reg=omit_reg,
                           only_first=only_first, n_best=n_best,
                           max_score=max_score, min_score=min_score,
                           only_chr=only_chr)
    elif file_format == 'score_tsv':
        regions = read_score_tsv(fp, omit_chr=omit_chr, omit_reg=omit_reg,
                                 only_first=only_first, n_best=n_best,
                                 max_score=max_score, min_score=min_score,
                                 only_chr=only_chr)
    else:
        raise ValueError(f'{file_format} is not a supported file format')
    
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


def read_bed(fp, omit_chr=(), omit_reg=(), only_first=False,
             only_chr=None, n_best=None, max_score=None, min_score=None):
    """Read and parse a bed file.

    Parameters
    ----------
    fp : str
        Path of the bed file.
    omit_chr : Iterable[str]
        An iterable containing names of chromosomes to be omitted.
    omit_reg : Iterable[str]
        An iterable containing names of regions to be omitted. This
        argument will be ignored if the bed file does not contain region
        names.
    only_first : bool
        Whether to keep only the first region when the bed file contains
        multiple regions with the same name. This argument will be
        ignored if the bed file does not contain region names.
    only_chr : Iterable[str]
        An iterable of chromosome names that will be exclusively
        considered.
    n_best : int
        Number of regions with the highest score to be returned.
    max_score : float
        Upper limit (inclusive) of scores. Only regions with scores
        lower or equal will be returned.
    min_score : float
        Lower limit (inclusive) of scores. Only regions with scores
        higher or equal will be returned.

    Returns
    -------
    List[Region]
        List of regions in the bed file
    """
    result = []
    omit_chr = set(omit_chr)
    omit_reg = set(omit_reg)
    seen_names = set()
    
    log_missing_info = False
    with open(fp) as f:
        for line in f:
            line = line.split()
            chromosome, start, end = line[:3]
            start = int(start)
            end = int(end)
            name = None
            score = 0
            pos_strand = True
            # positive strand is assumed if the bed file does not
            # contain the information
            
            try:
                name = line[3]
                score = float(line[4])
                pos_strand = line[5] == '+'
            except IndexError:
                log_missing_info = True
            
            # Omitting chromosomes
            if chromosome in omit_chr:
                continue

            # selecting chromosomes
            if only_chr and chromosome not in only_chr:
                continue

            # Omitting regions
            if name in omit_reg:
                continue
            
            # Keeping only first region
            if only_first:
                if name in seen_names:
                    continue
                elif name is not None:
                    seen_names.add(name)
            
            # flipping negative strands (does it make sense?)
            if pos_strand:
                region = Region(chromosome, start, end, name, pos_strand,
                                score)
            else:
                region = Region(chromosome, end, start, name, pos_strand,
                                score)
            
            result.append(region)

            # n highest score filter
            if n_best is not None:
                if n_best < len(result):
                    n_best_regions = sorted(result,
                                            key=lambda t: t.score,
                                            reverse=True)
                    n_best_regions = set(n_best_regions[:n_best])
                    result = list(
                        filter(lambda x: x in n_best_regions, result))

            # score range filter
            if max_score is not None:
                result = list(filter(lambda x: x.score <= max_score))
            if min_score is not None:
                result = list(filter(lambda x: x.score >= min_score))
                
    if log_missing_info:
        logging.info('Some optional fields were missing in the file')
    
    return result


def read_score_tsv(fp, omit_chr=(), omit_reg=(), only_first=None,
                   only_chr=None, n_best=None, max_score=None, min_score=None):
    """Read and parse a custom score tsv file.

        Parameters
        ----------
        fp : str
            Path of the bed file.
        omit_chr : Iterable[str]
            An iterable containing names of chromosomes to be omitted.
        omit_reg : Iterable[str]
            An iterable containing names of regions to be omitted. This
            argument will be ignored if the bed file does not contain
            region names.
        only_first : bool
            Whether to keep only the first region when the bed file
            contains multiple regions with the same name. This argument
            will be ignored if the bed file does not contain region
            names.
        only_chr : Iterable[str]
            An iterable of chromosome names that will be exclusively
             considered.
        n_best : int
            Number of regions with the highest score to be returned.
        max_score : float
            Upper limit (inclusive) of scores. Only regions with scores
            lower or equal will be returned.
        min_score : float
            Lower limit (inclusive) of scores. Only regions with scores
            higher or equal will be returned.

        Returns
        -------
        List[Region]
            List of regions in the file
        """
    result = []
    omit_chr = set(omit_chr)
    omit_reg = set(omit_reg)
    seen_names = set()
    
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
            
            # omitting chromosomes
            if chromosome in omit_chr:
                continue
            
            # omitting regions
            if name in omit_reg:
                continue
            
            # selecting chromosomes
            if only_chr and chromosome not in only_chr:
                continue
            
            # Keeping only first region
            if only_first:
                if name in seen_names:
                    continue
                elif name is not None:
                    seen_names.add(name)
            
            result.append(Region(chromosome, start, end, name, True, score))
    
    # n highest score filter
    if n_best is not None:
        if n_best < len(result):
            n_best_regions = sorted(result,
                                    key=lambda t: t.score,
                                    reverse=True)
            n_best_regions = set(n_best_regions[:n_best])
            result = list(filter(lambda x: x in n_best_regions, result))
    
    # score range filter
    if max_score is not None:
        result = list(filter(lambda x: x.score <= max_score))
    if min_score is not None:
        result = list(filter(lambda x: x.score >= min_score))
    
    return result
