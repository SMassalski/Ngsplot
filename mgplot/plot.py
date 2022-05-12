import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, LogNorm


# PROPOSED FEATURE: Joint plot
# DESIGN: Remove Axis parameters out of the function?
def average_profile(x, smooth=None, fig_kw=None, **plot_kw):
    r"""Plot the average profile of x
    
    Parameters
    ----------
    x : arr
        2D array containing the region signals.
    smooth : float, bool or None
        Smoothing factor to be used. If `True` std(profile)/2 will
        be used. If None (default) the plot won't be smoothed.
    fig_kw : dict
        Keyword arguments to be passed to plt.Figure.
    \*\*plot_kw : dict
        All additional keyword arguments will be passed to ax.set()

    Returns
    -------
    plt.Figure, axes.Axes
        Figure and axes of the plot.
    """
    if fig_kw is None:
        fig_kw = {}
        
    profile = np.mean(x, axis=0)
    
    # Smoothing
    if smooth is not None:
        if smooth is True:
            smooth = np.std(profile)/2
        space = np.arange(x.shape[1])
        spl = UnivariateSpline(space, profile, s=smooth)
        profile = spl(space)
        
    fig, ax = plt.subplots(**fig_kw)
    ax.plot(profile)
    ax.set(**plot_kw)
    
    return fig, ax


# TODO: Implement downsample
# FIXME: Somtimes overflow occurs
#   Right-most tick is not visible
def heatmap(x, sort=False, downsample=False, colorbar=False, fig_kw=None,
            imshow_kw=None, **plot_kw):
    r"""Plot the average profile of x
    
    Parameters
    ----------
    x : arr
        2D array containing the region signals.
    sort : bool
        Whether to sort the regions by sum(x[i]).
    downsample : int or None
        <<< Not implemented >>
        The factor by which the signals should be downsampled. If None
        (default) no downsampling will occur.
    colorbar : bool
        Whether to show a colorbar.
    fig_kw : dict
        Keyword arguments to be passed to plt.Figure.
    imshow_kw : dict
        Keyword arguments to be passed to imshow.
    \*\*plot_kw : dict
        All additional keyword arguments will be passed to ax.set()

    Returns
    -------
    plt.Figure, axes.Axes
        Figure and axes of the plot.
    """
    if fig_kw is None:
        fig_kw = {}
    if imshow_kw is None:
        imshow_kw = {}
        
    # Sorting
    if sort:
        idx = np.argsort(np.sum(x, axis=1))
        idx = np.flip(idx)
        x = np.take(x, idx, axis=0)
        
    if downsample:
        pass
    
    fig, ax = plt.subplots(**fig_kw)
    ax.set(**plot_kw)
    im = ax.imshow(x, **imshow_kw)
    
    if colorbar:
        fig.colorbar(im)
    
    return fig, ax


def make_ticks(roi, flank, body=None, n_ticks=2, roi_labels=None):
    """Construct ticks and tick labels for the plots
    
    Parameters
    ----------
    roi : str, one of ['start', 'end', 'body']
        Region of interest that will be included in the plot.
    flank : int
        Length of the flanking segments.
    body : int or None
        Length of the body segment. Only needed if `region_part` is set
        to 'body'.
    n_ticks : int
        Number of ticks per segment (left flank, right flank and body)
        in addition to start and end ticks.
    roi_labels : dict
        Mapping of labels to be used instead of 'start' and 'end'

    Returns
    -------
    ticks : list
        Positions of ticks.
    tick_labels : list
        Labels of the ticks.
        
    Raises
    ------
    ValueError
        If `region_part` == 'body' and `body` is None.
    """
    if roi_labels is None:
        roi_labels = {
            'start': 'start',
            'end': 'end'
        }
    roi_label = ''
    if roi != 'body':
        roi_label = 'TSS' if roi == 'start' else 'TSE'
    if body is None and roi == 'body':
        raise ValueError("Integer value for `body` is required when "
                         "`roi` == 'body'")
        
    if n_ticks == 0:
        ticks = [flank, flank + body] if roi == 'body' else [flank]
        tick_labels = ['TSS', 'TSE'] if roi == 'body' else [roi_label]
        
    elif roi == 'body':
        l_flank_labels = [str(x) for x in range(-flank, 0, flank // n_ticks)]
        l_flank_ticks = list(range(0, flank, flank // n_ticks))

        r_flank_labels = \
            list(range(flank // n_ticks, flank + 1, flank // n_ticks))
        r_flank_ticks = [flank + body + x for x in r_flank_labels]
        r_flank_labels = [f'+{x}' for x in r_flank_labels]

        body_labels = np.linspace(0, 1, n_ticks + 1, endpoint=False)[1:]
        body_labels = [f'{x :.2f}' for x in body_labels]
        body_ticks = np.linspace(flank, flank+body, n_ticks+1, endpoint=False)
        body_ticks = list(body_ticks[1:])
        
        ticks = l_flank_ticks \
            + [flank] \
            + body_ticks \
            + [flank + body] \
            + r_flank_ticks
        tick_labels = l_flank_labels \
            + [roi_labels['start']] \
            + body_labels \
            + [roi_labels['end']] \
            + r_flank_labels

    # Ticks for 'start' and 'end'
    else:
        l_flank_labels = list(range(-flank, 0, flank // n_ticks))
        r_flank_labels = \
            list(range(flank // n_ticks, flank + 1, flank // n_ticks))
        r_flank_labels = [f'+{x}' for x in r_flank_labels]
        tick_labels = l_flank_labels \
            + [roi_labels[roi]] \
            + r_flank_labels
        
        ticks = list(range(0, flank * 2 + 1, flank // n_ticks))
            
    return ticks, tick_labels


def make_norm(x, norm_type='lin', drop=1e-3):
    """Create norm for normalization of data to 0-1 range
    
    Parameters
    ----------
    x : arr
        Array containing the values that will be scaled.
    norm_type : str, one of ['lin', 'log']
        A string indicating the type of scale to use:
        lin: Linear scale (default)
        log: Logarithmic scale
    drop : float between 0 (inclusive) and 0.5 (exclusive)
         Value indicating quantiles that will be ignored when
         determining the scale.
         
         Data values < quantile(`drop`) and > quantile(1-`drop`) will be
         clipped.

    Returns
    -------
    Normalize object
        The norm.
    """
    if norm_type == 'lin':
        v_min = np.quantile(x, drop)
        v_max = np.quantile(x, 1-drop)
        norm = Normalize(v_min, v_max, True)
        
    elif norm_type == 'log':
        v_min = np.nanquantile(np.where(x > 0, x, np.nan), drop)
        v_max = np.quantile(x, 1 - drop)
        norm = LogNorm(v_min, v_max, True)
    else:
        raise ValueError(f'Unrecognized norm_type {norm_type}')
        
    return norm
    