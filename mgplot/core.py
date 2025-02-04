import logging
from collections import namedtuple
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from .util import group_by_chromosome


Region = namedtuple('Region', ['chromosome', 'start', 'end', 'name',
                               'positive_strand', 'score'])
SignalSegment = namedtuple('Signal', ['chromosome', 'start', 'end', 'value'])


# TODO:
#  * Tqdm capabilities.
def query(regions, signal, roi='start', flank=1000, body=None,
          flip_negative=False):
    """Construct signals for given regions and region of interest.

    Parameters
    ----------
    regions : list
        List of regions
    signal : list
        SignalSegments from which the chromosomes signal will be
        constructed.
    roi : str, one of ['start', 'body', 'end']
        Region of interest for which a signal will be returned.
        If 'body' is selected the regions will be normalized to a common
         length equal to `body` (final signal for each region will have
         a length of `body` + 2 * `flank`).
    flank : int
        Number of bp to be added to the start and end of each region.
    body : int or None
        Length to which the regions will be normalized to if selected
        region_part is 'body'. If None (default) body will be set to
        max(flank * 2, 1000)
    flip_negative : bool
        Whether to flip the signal of regions on the negative sense
        strand.

    Returns
    -------
    arr
        2D array where each row is the signal for a region from
        `regions`.

    Raises
    ------
    ValueError
        If `roi` is not one of the correct values.
    """
    signal = group_by_chromosome(signal)
    regions = group_by_chromosome(regions)
    
    # argument defaults
    if body is None and roi == 'body':
        body = max(flank * 2, 1000)
    
    # result array setup
    region_count = 0
    for regs in regions.values():
        region_count += len(regs)
    normalized_len = 2 * flank + body if roi == 'body' else 2 * flank + 1
    result = np.zeros((region_count, normalized_len), dtype=np.float16)
    
    idx = 0
    for chromosome, regs in regions.items():
        if chromosome not in signal:
            logging.warning(f'Signal info missing for chromosome {chromosome}')
            continue
        
        # create chromosome signal array
        
        # expanding chromosome signal by flanks to handle regions that after
        # expansion would reach outside the signal otherwise
        chr_signal = np.zeros(max(signal[chromosome], key=lambda x: x.end).end
                              + 2 * flank)
        for segment in signal[chromosome]:
            chr_signal[segment.start+flank:segment.end+flank] = segment.value
        
        for region in regs:
            
            # Start and end points for region parts
            if roi == 'start':
                l_flank_start = region.start
                l_flank_end = l_flank_start + flank
                r_flank_start = l_flank_end + 1
                r_flank_end = r_flank_start + flank
            elif roi == 'body':
                l_flank_start = region.start
                l_flank_end = l_flank_start + flank
                r_flank_start = region.end + flank + 1
                r_flank_end = r_flank_start + flank
            elif roi == 'end':
                l_flank_start = region.end
                l_flank_end = l_flank_start + flank
                r_flank_start = l_flank_end + 1
                r_flank_end = r_flank_start + flank
            else:
                raise ValueError("Incorrect `roi` value. Must be one"
                                 "of: ['start', 'body', 'end']")
            
            # Extracting data
            reg_signal = np.zeros(normalized_len)
            if flank > 0:  # Erroneous behavior for right flank when flank == 0
                reg_signal[:flank] = chr_signal[l_flank_start:l_flank_end]
                reg_signal[-flank:] = chr_signal[r_flank_start:r_flank_end]
            # DESIGN : is there a better way to plot body?
            #   Try out normalization of whole sequence
            if roi == 'body':
                reg_signal[flank:-flank or body] = \
                    normalize(chr_signal[l_flank_end:r_flank_start], body)
            else:
                reg_signal[flank] = chr_signal[l_flank_end]
            
            # Flipping regions on negative strand
            if flip_negative and not region.positive_strand:
                reg_signal = np.flip(reg_signal)
            
            result[idx] = reg_signal
            idx += 1
    return result


def normalize(arr, size):
    """Normalize an array to a specific length"""
    spl = InterpolatedUnivariateSpline(
        np.linspace(0, len(arr), len(arr)), arr, k=3)
    return spl(np.linspace(0, len(arr), size))
