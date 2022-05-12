import argparse
import sys
import logging
import numpy as np
from matplotlib.colors import Normalize, LogNorm
import seaborn as sb
from mgplot.io import load_regions, read_bedgraph
from mgplot.core import query
from mgplot.plot import heatmap, average_profile, make_ticks


# DESIGN: Keep defaults out of `add_argument` for argument hierarchy
# TODO:
#  * Update README
#  * Rename `region_part` to `roi`
#  * Add option for showing the plots
parser = argparse.ArgumentParser()
# Input files
parser.add_argument('-c', '--config_file', type=str,
                    help="Path of the file containing configuration arguments"
                         "and values")
parser.add_argument('-i', '--signal_file', type=str,
                    help="Path of the bedgraph file containing the signal")
parser.add_argument('-r', '--region_file', type=str,
                    help="Path of the file containing genomic regions")
parser.add_argument('-rf', '--reg_file_format', type=str, default='bed',
                    help="Format of the file containing regions "
                         "['bed','score_tsv']")
# DESIGN: Combine with signal_file (as input file)
parser.add_argument('-rp', '--replot', type=str,
                    help="File containing matrix to be replotted")

# DESIGN: Remove hm_file, avg_file, mat_file options?
# TODO: Implement out_prefix and out_dir
# Output files
parser.add_argument('-op', '--out_prefix', type=str,
                    help="Prefix to be used for output files")
parser.add_argument('-od', '--out_dir', type=str,
                    help="Path of the output directory")
parser.add_argument('-oh', '--hm_file', type=str,
                    help="Heatmap output filename")
parser.add_argument('-oa', '--avg_file', type=str,
                    help="Average profile output filename")
parser.add_argument('-om', '--mat_file', type=str,
                    help="Matrix output filename")

# DESIGN: Change to body, end, start. Plot should also use 'TSS', 'TSE'
#   and 'gene_body' labels
parser.add_argument('-p', '--region_part', type=str, default='TSS',
                    help="Region part to be plotted ['start','end','body']")
parser.add_argument('-fl', '--flank', type=int, default=3000,
                    help="Length of flanking fragments to be plotted with the"
                         " selected region part")
parser.add_argument('--body', type=int,
                    help="Length to which regions will be normalized if"
                         " region_part 'body' was selected")

# Filter options
# Elementwise filters
parser.add_argument('--only_chr', nargs='+', type=str,
                    help="The names of chromosomes to be exclusively"
                         " considered")
parser.add_argument('--omit_chr', nargs='+', type=str,
                    help="The names of chromosomes to be omitted")
parser.add_argument('--omit_reg', nargs='+', type=str,
                    help="Names of features to be ignored")
parser.add_argument('--min_score', type=float,
                    help="Ignore regions with scores lower than the given"
                         " value")
parser.add_argument('--max_score', type=float,
                    help="Ignore regions with scores higher than the given"
                         " value")
# Aggregate filters
parser.add_argument('--n_best', type=int,
                    help="The number of regions with the highest score"
                         " to be used")
parser.add_argument('--only_first', action='store_true',
                    help="Use only the first feature with"
                         " the same name")

# Plot options
# TODO: n_ticks help
parser.add_argument('--n_ticks', type=str,
                    help="Number of ticks per plot segment (left flank, right "
                         "flank and body) in addition to TSS and TSE ticks")
parser.add_argument('--plot_type', type=str,
                    help="Type of plot to be generated "
                         "['avg_prof','heatmap','both']; default = avg_prof")
parser.add_argument('--sort', action='store_true',
                    help="Whether to sort the matrix used for generating "
                         "a heatmap")
parser.add_argument('--cmap', type=str,
                    help="Colormap used in the heatmap; default = 'Reds'")
parser.add_argument('--smooth', type=str,
                    help="Smoothing factor used when smoothing the"
                         " average profile with a spline. If set to 'true' "
                         " <standard deviation> / 2 will be used. "
                         "Set to 0 for no smoothing.")
parser.add_argument('--hm_title', type=str, help="Title of the heatmap")
parser.add_argument('--avg_title', type=str,
                    help="Title of the average profile")
parser.add_argument('--cbar', action='store_true',
                    help="Show a color bar next to the heatmap")
parser.add_argument('--h_norm', type=str,
                    help="Type of norm to be used for the heatmap"
                         " colorscale ['lin','log']; default = 'lin'")
args = parser.parse_args()

config = {
            'region_part': 'TSS',
            'flank': '1000',
            'body': None,
            
            # Input options
            'signal_file': None,
            'region_file': None,
            'reg_file_format': 'bed',
            'replot': None,
            
            # Output options
            'hm_file': None,
            'avg_file': None,
            'mat_file': None,
            
            # Plot options
            'plot_type': 'avg',
            'cmap': 'Reds',
            'n_ticks': '1',
            'sort': False,
            'smooth': False,
            'hm_title': None,
            'avg_title': None,
            'h_norm': 'lin',
            'cbar': False,
            
            # Elementwise filters
            'omit_chr': None,
            'omit_reg': None,
            'only_chr': None,
            'min_score': None,
            'max_score': None,
            
            # Group filters
            'n_best': None,
            'only_first': False
        }


def read_config():
    if args.configFile:
        with open(args.configFile) as f:
            for i, line in enumerate(f):
                try:
                    k, v = line.split()
                except ValueError:
                    logging.warning(f'Wrong number of columns in line {i+1}'
                                    f'of the config file. Skipping...')
                    continue
                if k in config:
                    config[k] = v
    
    # override the config file values with the console values
    for arg in vars(args):
        if getattr(args, arg):
            config[arg] = getattr(args, arg)
    
    for arg in ['sort', 'only_first', 'cbar']:
        if type(config[arg]) == str:
            config[arg] = config[arg].lower() == 'true'
    for arg in ['omit_reg', 'omit_chr', 'only_chr']:
        if type(config[arg]) == str:
            config[arg] = config[arg].split(',')
    for arg in ['flank', 'n_ticks', 'n_best']:
        if type(config[arg]) == str:
            config[arg] = int(config[arg])
    for arg in ['min_score', 'max_score']:
        if type(config[arg]) == str:
            config[arg] = float(config[arg])
    for arg in ['region_part', 'plot_type', 'h_norm', 'reg_file_format']:
        if type(config[arg]) == str:
            config[arg] = config[arg].lower()


read_config()

if config['replot']:
    data = np.load(config['replot'])
elif config['region_file'] and config['signal_file']:
    regions = load_regions(config['region_file'],
                           file_format=config['reg_file_format'],
                           region_part=config['region_part'],
                           omit_chr=config['omit_chr'],
                           omit_reg=config['omit_reg'],
                           only_chr=config['only_chr'],
                           only_first=config['only_first'],
                           n_best=config['n_best'],
                           min_score=config['min_score'],
                           max_score=config['max_score'])
    signal = read_bedgraph(config['signal_file'])
    data = query(regions, signal, region_part=config['region_part'],
                 flank=config['flank'], body=config['body'])
    if config['mat_file']:
        np.save(config['mat_file'], data)
else:
    print("Input files not provided")
    sys.exit()
        
ticks, tick_labels = make_ticks(config['region_part'], config['flank'])
plot_kw = dict(xticks=ticks, xticklabels=tick_labels)
if config['plot_type'] in ['avg_prof', 'both']:
    with sb.axes_style("darkgrid"):
        ap_fig, ap_ax = average_profile(data, smooth=config['smooth'],
                                        title=config['avg_title'],
                                        **plot_kw)
    if config['avg_file'] is not None:
        ap_fig.save(config['avg_file'])
        
if config['plot_type'] in ['heatmap', 'both']:
    # TODO: Make norms work better.
    if config['h_norm'] == 'log':
        norm = LogNorm()
    else:
        norm = Normalize(data.min(), data.max())
    imshow_kw = dict(norm=norm, cmap=config['cmap'])
    hm_fig, hm_ax = heatmap(data, sort=config['sort'],
                            colorbar=config['colorbar'],
                            title=config['hm_title'],
                            imshow_kw=imshow_kw, **plot_kw)
    
    if config['hm_file'] is not None:
        hm_fig.save(config['hm_file'])
