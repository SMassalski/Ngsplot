import argparse
import sys
import logging
import numpy as np
import seaborn as sb
from mgplot.io import read_bedgraph, RegionParserBase, BedParser
from mgplot.core import query, Region
from mgplot.plot import heatmap, average_profile, make_ticks, make_norm
from mgplot.util import filter_regions


# DESIGN: Keep defaults out of `add_argument` for argument hierarchy
# TODO:
#  * Update README
#  * Add option for showing the plots
def get_cli_args(argv):
    parser = argparse.ArgumentParser()
    # Input files
    parser.add_argument('-c', '--config_file', type=str,
                        help="Path of the file containing configuration"
                             " arguments and values")
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
    
    parser.add_argument('--roi', type=str, default='start',
                        help="Region of interest to be plotted "
                             "['start','end','body']")
    parser.add_argument('-fl', '--flank', type=int, default=3000,
                        help="Length of flanking fragments to be plotted with"
                             " the selected region of interest")
    parser.add_argument('--body', type=int,
                        help="Length to which regions will be normalized if"
                             " roi='body' was selected")
    
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
    parser.add_argument('--n_ticks', type=int,
                        help="Number of ticks per plot segment (left flank,"
                             " right  flank and body) in addition to TSS and "
                             "TSE ticks")
    parser.add_argument('--plot_type', type=str,
                        help="Type of plot to be generated "
                             "['avg_prof','heatmap','both']; default = "
                             "avg_prof")
    parser.add_argument('--sort', action='store_true',
                        help="Whether to sort the matrix used for generating "
                             "a heatmap")
    parser.add_argument('--cmap', type=str,
                        help="Colormap used in the heatmap; default = 'Reds'")
    parser.add_argument('--smooth', type=str,
                        help="Smoothing factor used when smoothing the"
                             " average profile with a spline. If set to 'true'"
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
    return parser.parse_args(argv)


def main():
    args = get_cli_args(sys.argv[1:])
    c = Config()
    if args.config_file is not None:
        c.from_file(args.config_file)
    c.update(vars(args))

    if c.replot:
        data = np.load(c.replot)
    elif c.region_file and c.signal_file:
        if c.reg_file_format == 'bed':
            parser = BedParser()
        elif c.reg_file_format == 'score_tsv':
            parser = ScoreTSVParser()
            
        else:
            raise ValueError(f'Unrecognized file_format {c.reg_file_format}')

        regions = parser.parse(c.region_file)
        regions = filter_regions(regions, omit_chr=c.omit_chr,
                                 omit_reg=c.omit_reg, only_chr=c.only_chr,
                                 only_first=c.only_first, n_best=c.n_best,
                                 min_score=c.min_score, max_score=c.max_score)
        signal = read_bedgraph(c.signal_file)
        data = query(regions, signal, roi=c.roi,
                     flank=c.flank, body=c.body)
        if c.mat_file:
            np.save(c.mat_file, data)
    else:
        print("Input files not provided")
        sys.exit()

    ticks, tick_labels = make_ticks(c.roi, c.flank,
                                    roi_labels={'start': 'TSS', 'end': 'TSE'})
    plot_kw = dict(xticks=ticks, xticklabels=tick_labels)
    if c.plot_type in ['avg_prof', 'both']:
        with sb.axes_style("darkgrid"):
            ap_fig, ap_ax = average_profile(data, smooth=c.smooth,
                                            title=c.avg_title,
                                            **plot_kw)
        if c.avg_file is not None:
            ap_fig.save(c.avg_file)

    if c.plot_type in ['heatmap', 'both']:
        imshow_kw = dict(norm=make_norm(data, c.h_norm),
                         cmap=c.cmap)
        hm_fig, hm_ax = heatmap(data, sort=c.sort, colorbar=c.colorbar,
                                title=c.hm_title, imshow_kw=imshow_kw,
                                **plot_kw)
    
        if c.hm_file is not None:
            hm_fig.save(c.hm_file)


class Config:
    
    def __init__(self):
        # general
        self.roi = 'start'
        self.flank = 1000
        self.body = None
        
        # input
        self.signal_file = None
        self.region_file = None
        self.reg_file_format = 'bed'
        self.replot = None
        
        # output
        self.mat_file = None
        self.avg_file = None
        self.hm_file = None
        
        # filters
        self.omit_chr = None
        self.omit_reg = None
        self.only_chr = None
        self.min_score = None
        self.max_score = None
        self.n_best = None
        self.only_first = None
        
        # plotting general
        self.plot_type = 'avg'
        self.n_ticks = 1
        
        # average profile
        self.avg_title = None
        self.smooth = False
        
        # heatmap
        self.sort = False
        self.hm_title = None
        self.h_norm = 'lin'
        self.cmap = 'Reds'
        self.colorbar = False
    
    @staticmethod
    def convert_type(value):
        """Convert string value to implied type"""
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        
        if value.lower() == 'true':
            return True
        elif value.lower() == 'false':
            return False
        elif value.lower() == 'none':
            return None
        elif ',' in value:
            return value.split(',')
        else:
            return value
        
    def from_file(self, fp):
        """Set option values from file"""
        with open(fp) as f:
            for i, line in enumerate(f):
                try:
                    k, v = line.split()
                except ValueError:
                    logging.warning(f'Wrong number of columns in line {i+1}'
                                    f'of the config file. Skipping...')
                    continue
                if k in self.__dict__:
                    setattr(self, k, self.convert_type(v))
                
    def update(self, args):
        """Set option values from dict
        
        Parameters
        ----------
        args : dict
            Mapping of option name to value.
        """
        for arg, val in args:
            if val is None:
                continue
            if arg in ['smooth']:
                setattr(self, arg, self.convert_type(val))
            if type(arg) == bool:
                setattr(self, arg, val or getattr(self, arg))


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


if __name__ == '__main__':
    main()
