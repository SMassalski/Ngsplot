import argparse
from sys import stdout
import logging
from collections import namedtuple
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sb


pars = argparse.ArgumentParser()
pars.add_argument('-cf', '--configFile', type=str,
                  help="File containing configuration arguments and values")
pars.add_argument('-g', '--gfile', type=str, help="File containing genes")
pars.add_argument('-i', '--infile', type=str,
                  help="Bedgraph file containing the signal")
pars.add_argument('-gt', '--gfiletype', type=str,
                  help="Type of the file containing genes ['bed','score_tsv'];"
                       " default = 'bed'")
pars.add_argument('-rp', '--replot', type=str,
                  help="File containing matrix to be replotted")

pars.add_argument('-oh', '--hmfile', type=str, help="Heatmap output filename")
pars.add_argument('-oa', '--avgfile', type=str,
                  help="Average profile output filename")
pars.add_argument('-om', '--matfile', type=str, help="matrix output filename")

pars.add_argument('-r', '--region', type=str,
                  help="Region to be plotted ['TSS','TSE','genebody']; "
                       "default = 'TSS'")
pars.add_argument('-fl', '--flank', type=str,
                  help="Length of flanking fragments to be plotted with the"
                       " selected region; default = 3000")

pars.add_argument('-only', '--chrmonly', nargs='+', type=str,
                  help="The exact names of chromosomes to be"
                       " exclusively considered")
pars.add_argument('-co', '--chrmomit', nargs='+', type=str,
                  help="The exact names of chromosomes to be"
                       " exclusively considered")
pars.add_argument('-go', '--gomit', nargs='+', type=str,
                  help="Names of features to be ignored")
pars.add_argument('-nb', '--nbest', type=str,
                  help="The number of features with the best score to be used")
pars.add_argument('-sr', '--scorerange', nargs=2, type=str,
                  help="The score range from which the features are selected")
pars.add_argument('-of', '--ofirst', action='store_true',
                  help="Whether to use only the first feature with"
                       " the same name [no value]; default = False")

pars.add_argument('-nt', '--nticks', type=str)
pars.add_argument('-p', '--plottype', type=str,
                  help="Type of plot to be generated "
                       "['avgprof','heatmap','both']; default = avgprof")
pars.add_argument('-s', '--sort', action='store_true',
                  help="Whether to sort the matrix used for generating "
                       "a heatmap [no value]; default = False")
pars.add_argument('-cm', '--cmap', type=str,
                  help="Colormap used in the heatmap; default = 'Reds'")
pars.add_argument('-sm', '--smooth', type=str,
                  help="Smoothing factor used when smoothing the average"
                       " profile with a spline. Set to 'false' or 0 if you"
                       " don't want to smooth; default = flank*1e-4")
pars.add_argument('-ht', '--hmtitle', type=str, help="Title of the heatmap")
pars.add_argument('-at', '--avgtitle', type=str,
                  help="Title of the average profile")
pars.add_argument('-cb', '--cbar', action='store_true',
                  help='Whether to show a colorbar next to the heatmap'
                       ' [no value]; default = False')
pars.add_argument('-hn', '--hnorm', type=str,
                  help="Type of norm to be used for the heatmap"
                       " colorscale ['lin','log']; default = 'lin'")
args = pars.parse_args()

config = {  # default argument values are stored here
    'region': 'TSS',
    'flank': '1000',
    
    'infile': None,
    'gfile': None,
    'gfiletype': 'bed',
    'replot': None,
    'dtype': 'int8',
    
    'hmfile': None,
    'avgfile': None,
    'matfile': None,
    
    'plottype': 'avgprof',
    'cmap': 'Reds',
    'nticks': '1',
    'sort': False,
    'smooth': False,
    'hmtitle': None,
    'avgtitle': None,
    'hnorm': 'lin',
    'cbar': False,
    
    'chrmomit': None,
    'chrmadd': None,
    'gomit': None,
    'chrmonly': None,
    
    'nbest': None,
    'scorerange': None,
    'ofirst': False
}


def read_config():
    if args.configFile:
        f = open(args.configFile, 'r')
        
        for line in f:
            l = line.split()
            if len(l) == 2:
                if l[0] in config:
                    config[l[0]] = l[1]
        f.close()
    
    # override the config file values with the console values
    for arg in vars(args):
        if getattr(args, arg):
            config[arg] = getattr(args, arg)
    
    for arg in ['sort', 'ofirst', 'cbar']:
        if type(config[arg]) == str:
            config[arg] = config[arg].lower() == 'true'
    for arg in ['gomit', 'chrmomit', 'gindx', 'scorerange', 'chrmonly']:
        if type(config[arg]) == str:
            config[arg] = config[arg].split(',')
    for arg in ['flank', 'nticks', 'nbest']:
        if type(config[arg]) == str:
            config[arg] = int(config[arg])
    for arg in ['region', 'plottype', 'hnorm', 'gfiletype']:
        if type(config[arg]) == str:
            config[arg] = config[arg].lower()


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


def query(fname, genome):
    f = open(fname, 'r')
    result = []
    line = f.readline().split()
    
    while line:
        chrm = line[0]
        if chrm in genome:
            
            signal = np.zeros(int(line[1]), dtype='int8')
            
            start = 0
            end = int(line[1])
            
            i = 1
            for gene in genome[chrm]:
                stdout.write(chrm + '\t' + str(i) + '\r')
                i += 1
                # Expansion of queried sequence by flanking sequences
                if region == 'genebody':
                    gStart = min(gene) - flank
                    gEnd = max(gene) + flank + 1
                else:
                    gStart = gene[TSE] - flank
                    gEnd = gene[TSE] + flank + 1
                
                if gStart >= end:
                    while line and int(line[1]) < gStart:
                        if line[0] != chrm:
                            break
                        x, y, z = int(line[1]), int(line[2]), int(line[3])
                        line = f.readline().split()
                    if gStart < y:
                        signal = np.append(
                            np.array([z] * (y - gStart), dtype='int8'),
                            np.zeros(int(line[1]) - y, dtype='int8'))
                    else:
                        signal = np.zeros(int(line[1]) - gStart, dtype='int8')
                    start = gStart
                    end = int(line[1])
                
                elif gStart > start:
                    signal = signal[gStart - start:]
                    start = gStart
                
                elif gStart < start and start == 0:
                    signal = np.append(np.zeros(start - gStart, dtype='int8'),
                                       signal)
                    start = gStart
                
                if gEnd + 1 > end:
                    signal = np.append(signal,
                                       np.zeros(gEnd - end, dtype='int8'))
                    x, y, z = int(line[1]), int(line[2]), int(line[3])
                    signal[x - start:y - start] = z
                    end = y
                    line = f.readline().split()
                    while gEnd + 1 > end:
                        if not line or line[0] != chrm:
                            end = gEnd
                            break
                        x, y, z = int(line[1]), int(line[2]), int(line[3])
                        signal[x - start:y - start] = z
                        end = y
                        line = f.readline().split()
                    if gEnd < x:
                        signal = np.append(signal,
                                           np.zeros(x - gEnd, dtype='int8'))
                    if gEnd < y:
                        signal = np.append(signal, np.array([z] * (y - x),
                                                            dtype='int8'))
                
                if gene[0] > gene[1]:
                    res = np.flip(signal[:gEnd - start], axis=0)
                else:
                    res = signal[:gEnd - start]
                
                if region == "genebody":
                    res = np.append(res[:flank],
                                    [normalize(res[flank:-flank], flank),
                                     res[-flank:]])
                
                result.append(res)
            
            print(chrm, 'done', i - 1)
            
            while line and line[0] == chrm:
                line = f.readline().split()
        else:
            while line and line[0] == chrm:
                line = f.readline().split()
    f.close()
    return np.array(result, dtype='int8')


def normalize(arr, size):
    spl = interpolate.InterpolatedUnivariateSpline(
        np.linspace(0, len(arr), len(arr)), arr, k=3)
    return spl(np.linspace(0, len(arr), size))


def plot(values):
    # plot vars set up
    genebody = config['region'] == 'genebody'
    plottype = config['plottype']
    
    if genebody:
        size = 3 * flank
        s = flank / 10000
    else:
        size = 2 * flank + 1
        s = flank / 10000
    
    nticks = config['nticks']
    if genebody:
        ticks = ['-' + str(i * flank // nticks) for i in
                 reversed(range(1, nticks + 1))] + ['TSS'] + [
                    '%.2f' % float(i / nticks) for i in range(1, nticks)] + [
                    'TSE'] + ['+' + str(i * flank // nticks) for i in
                              range(1, nticks + 1)]
        if nticks == 0:
            tickvals = [flank, 2 * flank]
        else:
            tickvals = [i * flank // nticks for i in range(0, nticks)] + [
                flank + i * flank // nticks for i in range(0, nticks)] + [
                           2 * flank + i * flank // nticks for i in
                           range(0, nticks + 1)]
    else:
        ticks = ['-' + str(i * flank // nticks) for i in
                 reversed(range(1, nticks + 1))] + [region.upper()] + [
                    '+' + str(i * flank // nticks) for i in
                    range(1, nticks + 1)]
        if nticks == 0:
            tickvals = [flank]
        else:
            tickvals = np.linspace(0, size, 2 * nticks + 1)
    
    if config['smooth']:
        if config['smooth'].lower() == 'true':
            smooth = flank / 10000
        elif config['smooth'].lower() == 'false':
            config['smooth'] = False
        else:
            smooth = float(config['smooth'])
    
    # plot avgprof
    if plottype in ['avgprof', 'both']:
        
        print('calculating mean...')
        
        if config['smooth']:
            spl = interpolate.UnivariateSpline(np.linspace(0, size, size),
                                               np.mean(values, axis=0,
                                                       dtype='float16'),
                                               s=smooth)
            avgprof = spl(np.linspace(0, size, size))
        
        else:
            avgprof = np.mean(values, axis=0, dtype='float16')
        
        print('plotting...')
        with sb.axes_style("darkgrid"):
            plt.plot(np.linspace(0, size, size), avgprof)
            
            plt.xticks(tickvals, ticks)
            plt.xlim((0, size))
            plt.ylim((avgprof.min() - 0.05, avgprof.max() * 1.05))
            if config['avgtitle']:
                plt.title(config['avgtitle'])
            
            if config['avgfile']:
                plt.savefig(config['avgfile'])
            else:
                plt.show()
    # plot heatmap
    if plottype in ['heatmap', 'both']:
        
        if config['sort']:
            print('sorting...')
            b = np.sum(values, axis=1) * -1
            indx = b.argsort()
            values = np.take(values, indx, axis=0)
        
        if plottype == 'both':
            plt.figure()
        
        print('plotting...')
        with sb.axes_style("ticks"):
            if config['hnorm'] == 'lin':
                plt.imshow(values, aspect='auto', cmap=config['cmap'],
                           norm=colors.Normalize(vmin=values.min(),
                                                 vmax=values.max() * .8))
            elif config['hnorm'] == 'log':
                plt.imshow(values, aspect='auto', cmap=config['cmap'],
                           norm=colors.LogNorm())
            
            plt.xticks(tickvals, ticks)
            plt.xlim((0, size))
            
            if config['cbar']:
                plt.colorbar()
            if config['hmtitle']:
                plt.title(config['hmtitle'])
            
            if config['hmfile']:
                plt.savefig(config['hmfile'])
            else:
                plt.show()


Region = namedtuple('Region', ['chromosome',
                               'start',
                               'end',
                               'name',
                               'positive_strand'])


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
            Upper limit (inclusive) of scores. Only regions with scores lower
            or equal will be returned.
        min_score : float
            Lower limit (inclusive) of scores. Only regions with scores higher
            or equal will be returned.

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


# DESIGN: Filter after all regions are collected?
#   Implement a parser class.
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
        Upper limit (inclusive) of scores. Only regions with scores lower
        or equal will be returned.
    min_score : float
        Lower limit (inclusive) of scores. Only regions with scores higher
        or equal will be returned.

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


read_config()
region = config['region']
flank = config['flank']
TSE = region == 'tse'

if config['replot']:
    print('loading file')
    plot(np.load(config['replot']))
elif config['gfile'] and config['infile']:
    s = load_regions(config['gfile'], file_format=config['gfiletype'],
                     region_part=config['region'], omit_chr=config['chrmomit'],
                     omit_reg=config['gomit'], only_chr=config['chrmonly'],
                     only_first=config['ofirst'], n_best=config['nbest'],
                     min_score=float(config['scorerange'][0]),
                     max_score=float(config['scorerange'][1]))
    t = query(config['infile'], s)
    if config['matfile']:
        np.save(config['matfile'], t)
    print('array generated')
    plot(t)
    print('Done')
else:
    print("Input files not provided")
