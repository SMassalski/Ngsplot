# MGPLOT

MGPLOT is a tool for visualizing genomic feature data.  The script generates average profiles and heatmaps of selected genomic features based on a provided bedgraph file.


## Dependencies ##
   ---
1. Python (tested on 3.8)
2. NumPy (tested on 1.22.3)
3. SciPy (tested on 1.8.0)
4. Matplotlib (tested on 3.5.1)
5. Seaborn (tested on 0.11.2)

# Usage #
  ---
1. [Input Options](https://github.com/SMassalski/mgplot#input-options)

2. [Plot Options](https://github.com/SMassalski/mgplot#plot-options)

3. [Region Declaration](https://github.com/SMassalski/mgplot#region-declaration)

4. [Output](https://github.com/SMassalski/mgplot#output)

5. [Detailed Feature Selection](https://github.com/SMassalski/mgplot#detailed-feature-selection)
## Input Options 
|        Argument        |                      Description                      |  Default   |          Values           |
|:----------------------:|:-----------------------------------------------------:|:----------:|:-------------------------:|
| `-r ` `--region_file ` |            File containing genomic regions            |    None    |     file name or path     |
|  `-i` `--signal_file`  |          Bedgraph file containing the signal          |    None    |     file name or path     |
| `-c ` `--config_file`  |                  Configuration file                   |    None    |     file name or path     |
|      ` --format`       |         Format of the file containing regions         |    bed     |      bed , score_tsv      |
|       `--replot`       | Numpy binary file containing a matrix to be replotted |    None    |     file name or path     |


**Console arguments override configuration file arguments**

For the configuration file to be read correctly it should have only two whitespace separated columns. In each row the first column is the argument, the second is the value.
When an argument can take multiple values such as in `--omit_chr` and `--omit_reg` the values should be in the second column separated by commas.

Example of a correctly formatted configuration file:
```
roi         body
plot_type   both
n_ticks     2
sort        true
only_chr    chr1,chr2,chrX
min_score   1000,
max_score   10000
```


`--replot` allows you to plot previously saved data (using `--save_data`) without reading other input files again, with different plot settings. When using replot the `roi`, `body` and `flank` arguments must remain unchanged.


## Plot Options
|    Argument    |                                                     Description                                                      | Default  |                         Values                         |
|:--------------:|:--------------------------------------------------------------------------------------------------------------------:|:--------:|:------------------------------------------------------:|
| ` --plot_type` |                                             Type of plot to be generated                                             | avg_prof |               avg_prof , heatmap , both                |
|  `--n_ticks`   | Number of ticks on the x-axis per plot segment (left flank, right flank and body) in addition to start and end ticks |    1     |                 a non negative integer                 |
|    `--cmap`    |                                    Colormap to be used when generating a heatmap                                     |   Reds   |                matplotlib colormap name                |
|    `--sort`    |                                    Sort the matrix used for generating a heatmap                                     |  False   | no value in console,  'true' or 'false' in config file |
|   `--smooth`   |         Smoothing factor for average profile smoothing. Set to True to automatically determine based on std          |  False   |            boolean or a non negative float             |
|    `--norm`    |        Type of normalization to use (linear, logarithmic or symmetric logarithmic) for the heatmap colorscale        |   lin    |                  lin , log , sym_log                   |
|  `--hm_title`  |                                                 Title of the heatmap                                                 |   None   |                         string                         |
| `--avg_title`  |                                             Title of the average profile                                             |   None   |                         string                         |
|  `--colorbar`  |                                         Show a colorbar next to the heatmap                                          |  False   | no value in console, 'true' or 'false' in config file  |


When `--smooth` is set to true the curve is smoothed with a spline with a smoothing parameter equal to `std(mean(input_array))/2`.

Example average profile with `--smooth` set to 0.08:
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/avgprof_example3a.png "smooth = 0.08")

The same average profile with `--smooth` set to 0:
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/avgprof_example3b.png "smooth = 0")

Matplotlib colormap names can be found [here.](https://matplotlib.org/stable/gallery/color/colormap_reference.html)

`--norm` allows you to use a logarithmic or symmetric logarithmic normalization when linear normalization generates a poorly visible heatmap.

Example heatmap with `--norm` set to 'lin':
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/heatmap_example2c.png "norm = lin")

The same heatmap with `--norm` set to 'sym_log':
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/heatmap_example2b.png "norm = log")

## Region Declaration
| Argument  |                    Description                     | Default |          Values          |
|:---------:|:--------------------------------------------------:|:-------:|:------------------------:|
|  `--roi`  |          Region of interest to be plotted          |  start  |    start , end , body    |
| `--flank` |           Length of the flanking regions           |  1000   |   non-negative integer   |
| `--body`  | Length to which genomic regions will be normalized |  None   | positive integer or none |

The script will generate plots of the selected location with surrounding flanking segments.

When `--roi` is set to body the regions will be normalized with a cubic spline to have length equal to the value of `--body`.
If `--body` is set to None `max(2*flank, 1000)` will be used

## Output
|    Argument    |               Description                |       Default       |                         Values                         |
|:--------------:|:----------------------------------------:|:-------------------:|:------------------------------------------------------:|
| `--out_prefix` |    Prefix to be used for output files    |       mgplot        |                          str                           |
|  `--out_dir`   |       Path of the output directory       | current working dir |                          path                          |
| `--out_format` |       Format of plot output files        |         png         |                     png, pdf, svg                      |
| `--save_data`  | Save data for replotting (in a npy file) |        False        | no value in console,  'true' or 'false' in config file |

 
The `--save_data` argument allows you to save a matrix for replotting data without repeating data extraction. The matrix is saved in a numpy binary file.
 
 ## Detailed Feature Selection
|    Argument    |                                 Description                                  | Default |                        Values                         |
|:--------------:|:----------------------------------------------------------------------------:|:-------:|:-----------------------------------------------------:|
|  `--omit_reg`  |                        Genomic regions to be omitted                         |  None   |                     region names                      |
|  `--omit_chr`  |                          Chromosomes to be omitted                           |  None   |                   chromosome names                    |
|  `--only_chr`  |                   Chromosomes to be considered exclusively                   |  None   |                   chromosome names                    |
| `--min_score`  | Lower limit (inclusive) of score values of which features will be considered |  None   |                         float                         |
| `--max_score`  | Upper limit (inclusive) of score values of which features will be considered |  None   |                         float                         |
| `--only_first` |         Use only the first occurrence of regions with the same name          |  False  | no value in console, 'true' or 'false' in config file |
|   `--n_best`   |           Number of regions with the best scores to be considered            |  None   |               positive integer or none                |


The `--n_best` argument selects genes after chromosomes and genes are excluded.

# Example
---
The configuration file (example1.cfg):
```
region_file RNAseq_counts.tsv
signal_file H3K4me3.bedgraph
format  score_tsv
flank   2000
region  start
n_ticks 2
plot_type   both
smooth  true
save_data   true
norm    log
sort    false
n_best  2000
```
Running the script will generate an average profile and heatmap of H3K4me3 enrichment around the transcription start site on the first 2000 most expressed genes in the tsv file. The matrix will be saved for replotting as mgplot_data.npy.
```
python make_plots.py -c example1.cfg
```
The result:

![alt_text](https://github.com/SMassalski/mgplot/blob/master/Examples/example1/example1b.png "average profile")

![alt_text](https://github.com/SMassalski/mgplot/blob/master/Examples/example1/example1c.png "heatmap")

Replotting the heatmap with `--sort` set to True:
```
python make_plots.py -c example1.cfg --replot mgplot_data.npy --sort
```

![alt_text](https://github.com/SMassalski/mgplot/blob/master/Examples/example1/example1a.png "sorted heatmap")
