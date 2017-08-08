# MGPLOT

MGPLOT is a script for visualizing metagenomic data.


## Dependencies ##
   ---
1. Python (tested on 3.6.1)
2. NumPy (tested on 1.12.1)
3. SciPy (tested on 0.19.0)
4. Seaborn (tested on 0.7.1)

# Usage #
  ---
1. [Input Options](https://github.com/SMassalski/mgplot#input-options)

2. [Plot Options](https://github.com/SMassalski/mgplot#plot-options)

3. [Region Declaration](https://github.com/SMassalski/mgplot#region-declaration)

4. [Output](https://github.com/SMassalski/mgplot#output)

5. [Detailed Feature Selection](https://github.com/SMassalski/mgplot#detailed-feature-selection)
## Input Options 
| Argument              | Description                          | Default  | Values   |
|:---------------------:|:------------------------------------:|:--------:|:--------:|
|`-g ` `--gfile `| bed or tsv file containing information about the genomic features  |   None   | file name or path |
|`-i` `--infile`| bedgraph file |   None   | file name or path |
|`-gt`<br>` --gfiletype`| type of the gfile  | 'bed' | 'bed' , 'scoretsv' |
|`-cf `<br>`--configFile`| configuration file |   None   | file name or path |
|`-gi ` `--gindx`| indicies representing the location of fields in a bed file |   0,1,2,3,-1   | sequence of five integers |
|`-rp ` `--replot`| numpy binary file containing a matrix to be replotted |   None   | file name or path |

**Input files cannot have headers.**

**Console arguments override configuration file arguments**

For the configuration file to be read correctly it should have only two columns separated by TAB. In each row the first column is the argument, the second is the value.
When an argument can take multiple values such as in `--chrmomit` and `-gomit` the values should be in the second column separated by commas.
Example of a correctly formatted configuration file:
```
region      genebody
plottype    both
nticks      2
sort        true
chrmonly    chr1,chr2,chrX
scorerange  1000,10000
```

The `--gindx` sequence indicates in which columns in the `.bed` file the script will find information about chromosome name, start position, end position, feature name, strand respectively. 

`--replot` allows you to plot previously saved data (using `--matfile`) without reading other input files again, with different plot settings. When using replot the _region_ and _flank_ arguments must remain unchanged.


## Plot Options
| Argument | Description| Default  | Values   |
|:----------:|:----------:|:--------:|:--------:|
|`-p` <br>` --plottype` | type of plot to be generated |   'avgprof'   |  'avgprof' , 'heatmap' , 'both'  |
|`-nt` <br>`--nticks`| number of ticks on the x-axis  |   1   | a non negative integer |
|`-cm` `--cmap`| colormap to be used when generating a heatmap   | 'Reds' | matplotlib colormap name |
|`-s` `--sort`| whether to sort the matrix used for generating a heatmap |   False   | no value in console,  'true' or 'false' in config file |
|`-sm`<br> `--smooth`| whether to smooth the average profile curve before ploting |   False  | boolean or a non negative float  |
|`-hn` `--hnorm`| whether to use linear or symmetric logarithmic normalization for the colorscale  | 'lin' | 'lin' , 'log' |
|`-ht` <br>`--hmtitle`| title of the heatmap  |   None  | string |
|`-at` <br>`--avgtitle`| title of the average profile  |   None  | string |
|`-cb` <br>`--cbar`| whether to show a colorbar next to the heatmap  |   False   | no value in console,<br> 'true' or 'false' in config file |


`--nticks` determines how many ticks surround the plotted region. For TSS and TSE *nticks* is the number of ticks on each side of the region of interest. For genebody it is the number of ticks left of TSS, between TSS and TSE (nticks-1), and right of TSE.
For example a run with argument `--nticks 2` will generate plots with 5 ticks for TSS and TSE, and 7 ticks for genebody.

When `--smooth` is set to true the curve is smoothed with a spline with a smoothing parameter equal to flank_length * 1e-4. To change the value of the smoothing parameter set the value of `--smooth` to a float.

Example average profile with `--smooth` set to 0.08:
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/avgprof_example3a.png "smooth = 0.08")

The same average profile with `--smooth` set to 0:
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/avgprof_example3b.png "smooth = 0")

Matplotlib colormap names can be found [here.](https://matplotlib.org/users/colormaps.html)

`--hnorm` allows you to use a symmetric logarithmic normalization when linear normalization generates a poorly visable heatmap. The 'log' heatmaps use gaussian interpolation

Example heatmap with `--hnorm` set to 'lin':
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/heatmap_example2c.png "hnorm = lin")

The same heatmap with `--hnorm` set to 'log':
![alt text](https://github.com/SMassalski/mgplot/blob/master/Examples/heatmap_example2b.png "hnorm = log")

## Region Declaration
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-r` `--region`| location of the genomic regions to be plotted|'TSS'| 'TSS' , 'TSE' , 'genebody'|
|`-fl` `--flank` |length of the flanking regions|1000| positive integer|

The script will generate plots of the selected location with surrounding flanking regions.

When `--region` is set to genebody the regions will be normalized with a cubic spline to have length equal to the value of `--flank`.

## Output
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-oa` `--avgfile`|average profile output file|None|file name or path|
|`-oh` `--hmfile`|heatmap output file|None|file name or path|
|`-om` `--matfile`|matrix output file|None|file name or path|

 The format of the output file will be deduced from the extension of the provided filename. The supported formats depend on the matplotlib backend you're using. Most backends support *png, pdf, ps, eps* and *svg*. If a output filename is provided the corresponding plot will not be displayed.
 
The `--matfile` argument allows you to save a matrix for replotting data without reading other input files again. The matrix is saved in a numpy binary file. A `.npy` extension will be appended if the provided filename does not have one.
 
 ## Detailed Feature Selection
| Argument | Description| Default  | Values   |
|:--------:|:----------:|:--------:|:--------:|
|`-co` <br> `--chrmomit`|chromosomes to be omitted|None|chromosome names|
|`-only`<br> `--chrmonly`|chromosomes to be considered exclusively|None|chromosome names|
|`--go` <br>`--gomit`|genomic regions to be omitted|None|first characters of feature names|
|`-nb`<br> `--nbest`|number of regions with the best scores to be considered (tsv files)|None|positive integer|
|`-sr` <br>`--scorerange`|range of score values from which features should be considered (tsv file)|None| two numbers seperated by  a comma|
|`-of`<br> `--ofirst`|whether to use only the first occurence of a region in a bed file|False|no value in console,<br> 'true' or 'false' in config file |

The `--gomit` argument makes the script omit all genomic features that **begin with** one of the given strings.

The `--nbest` argument selects genes after chromosomes and genes are excluded.

When using `--scorerange` the first value should be smaller than the other.

# Examples
---

