# filtering
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
    result = list(filter(lambda x: keep_region(x, **kwargs), regions))
    
    # only_first filter
    
    if only_first:
        seen_names = set()
        tmp_result = []
        for r in result:
            if r.name in seen_names:
                continue
            else:
                if r.name is not None:
                    seen_names.add(r.name)
                tmp_result.append(r)
        result = tmp_result
    
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


def keep_region(region, omit_chr=(), omit_reg=(), only_chr=None,
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
    if max_score is not None and region.score < max_score:
        return False
    if min_score is not None and region.score > min_score:
        return False
    
    return True


def group_by_chromosome(x):
    """Group a collection of SignalSegments or Regions by chromosome
    
    Parameters
    ----------
    x : Iterable
        Collection of SignalSegments or Regions.

    Returns
    -------
    dict
        Mapping of chromosome names to corresponding elements
    """
    result = {}
    for i in x:
        if i.chromosome not in result:
            result[i.chromosome] = []
        result[i.chromosome].append(i)
    return result
