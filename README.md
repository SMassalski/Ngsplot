# PMGPLOT

## Dependancies
1. Python 3.6.1
2. NumPy 1.12.1
3. SciPy 0.19.0

# Usage


## Input Options
| Argument              | Description                          | Default  | Values   |
|:---------------------:|:------------------------------------:|:--------:|:--------:|
|```-g --gfile     ```| bed or tsv file containing information about the genes  |   None   | file name or path |
|```-i --infile    ```| bedgraph file containing the signal  |   None   | file name or path |
|```-gt --gfiletype```| type of file containing information about the genes   | 'bed'    | [ 'bed' , 'scoretsv' ]|
|```-cf --configFile```| configuration file |   None   | file name or path |

All input files cannot have headers.

For the configuration file to be read correctly it should have only two columns separated by TAB. In each row the first column is the argument, the second is the value.
When an argument can take multiple values such as in `--chrmomit` and `-gomit` the values should be in the second column separated by commas.
Example of a configuration file:
```
region  genebody
plottype    both
nticks  2
sort    true
chrmomit chr1,chrX,chrY
scorerange 1000,10000
```


## Plot Options
| Argument | Description| Default  | Values   |
|:----------:|:----------:|:--------:|:--------:|
|`-p` ` --plottype` | type of plot to be generated |   'avgprof'   | [ 'avgprof' , 'heatmap' , 'both' ] |
|`-nt` `--nticks`| number of ticks on the x-axis  |   1   | a non negative integer |
|`-cm` `--cmap`| colormap to be used when generating a heatmap   | 'Reds' | matplotlib colormap name |
|`-s` `--sort`| whether to sort the matrix used for generating a heatmap |   False   | no value in console, [ 'true' , 'false' ] in config file |
|`-sm` `--smooth`| whether to smooth the average profile curve before ploting |   False  | boolean or a non negative float  |

`--nticks` determines how many ticks surround the plotted region. For TSS and TSE *nticks* is the number of ticks on each side of the region of interest. For genebody it is the number of ticks left of TSS, between TSS and TSE, and right of TSE.
For example a run with argument `--nticks 1` will generate plots with 3 ticks for TSS and TSE, and 5 ticks for genebody.

When `--smooth` is set to true the curve is smoothed with a splin with a smoothing parameter equal to flank_length * 1e-4. To change the value of the smoothing parameter set the value of `--smooth` to a float.

Matplotlib colormap names can be found [here.](https://matplotlib.org/users/colormaps.html)

## Region Declaration
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-r` `--region`| genomic region to be plotted|'TSS'|[ 'TSS' , 'TSE' , 'genebody' ]
|`-fl` `--flank` |length of the flanking regions|1000| positive integer|


## Output
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-oa` `--avgfile`|average profile output file|None|file name or path|
|`-oh` `--hmfile`|heatmap output file|None|file name or path|

 The format of the output file will be deduced from the extension of the provided filename. The supported formats depend on the matplotlib backend you're using. Most backends support *png, pdf, ps, eps* and *svg*. If a output filename is provided the corresponding plot will not be displayed.
 
 ## Detailed Gene Selection
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-co` <br> `--chrmomit`|chromosomes to be omitted|None|chromosome names|
|`-ca` `--chrmadd`|additional chromosomes to be considered|None|chromosome names|
|`--go` `--gomit`|genes to be omitted|None|first characters of gene names|
|`-nb` `--nbest`|number of genes with the best scores to be considered (tsv files)|None|positive integer|
|`-sr` `--scorerange`|range of score values from which genes should be considered (tsv file)|None| two numbers seperated by  a comma|
|`-of` `--ofirst`|whether to use only the first occurence of a gene in a bed file|False|no value in console,<br> [ 'true' , 'false' ] in config file |

The `--gomit` argument makes the script omit all genes that **begin with** one of the given strings.
The `--nbest` argument selects genes after chromosomes and genes are excluded.


