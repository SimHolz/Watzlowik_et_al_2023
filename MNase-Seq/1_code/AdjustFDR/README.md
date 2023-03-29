# AdjustFDR.R

The way danpos calculates FDR values is in my opinion entirely not correct! Therefore I wrote `AdjustFDR.R`, an R-script for correct determination of FDRs for nucleosome dynamics statistics. It simulates additional, random nucleosomes, calculates their pvalue and uses those simulated values for correct FDR estimation.

The number of simulated values for FDR calculation can be adjusted. If not provided, danpos detected nucleosomes x 100 pvalues will be simulated. 5000000 simulated pvalues take about 40-50min per sample. This scales nearly linearly with simulation n.

For estimation of fuzziness_diff_fdr & smt_diff_fdr you can choose between 3 types of simulation strategies:

- "same" samples a random position and calculates pvalue for score between treatment and control on this position. Used in original DANPOS.
- "random" samples two random positions and calculates pvalue for treatment score on one position and control score on the other position.
- "near" samples a random positions and a second position near (+-150bp) that position. It calculates pvalue for treatment score on one position and control score on the other position. (default)

If DANPOS wig files have not been converted to bigwigs yet, a chromosome.sizes file has to be provided and bigwig files will be generated.

## to-do:

The main functionality of the script is there but it is still under development. Please ask me if there are any questions or it doesn't work.

Currently there are the following caveats:

- [x] ~~you still have to change the path of danpos output manually in the script~~

- [x] ~~pooled wig files have to be converted to bigwig first.~~

- [ ] danpos wiggle step size has to be 1 (Danpos option `-a 1` )

- [x] ~~handling of NAs in control_smt_loca is not ideal.~~

This might be changed in future versions!

## Dependancies:
The script is dependant on the following R packages:

- `argparser`
- `tidyverse`
- `plyranges`
- `fs`

The packages can be installed via CRAN (`argparser`, `tidyverse` and `fs`) and bioconductor (`plyranges`).

## Command Line Usage:

    usage: 
    AdjustFDR.R [--help] [--output_path OUTPUT_PATH] [--chr_sizes CHR_SIZES] [--n_simulations N_SIMULATIONS] 
    [--sampling_type_occupancy SAMPLING_TYPE_OCCUPANCY] [--sampling_type_fuzziness SAMPLING_TYPE_FUZZINESS]
    [--sample_names SAMPLE_NAMES] [--control_names CONTROL_NAMES] input_path

### Arguments

#### input_path

Input path. Must be a path to DANPOS output.

#### -h, --help

show this help message and exit

#### -o, --output_dir [default: /adjustedFDR/]

Output directory.

#### -c, --chr_sizes

Chrom.sizes file for wig to bigwig conversion.

#### -n, --n_simulations [default: 100x fdr-values to estimate]

Number of simulated p-values for FDR estimation. 

#### -sto, --sampling_type_occupancy [default: near]

sampling type for fdr estimation of occupancy change fdr. 
Either 'same', 'random' or 'near'.

#### -stf, --sampling_type_fuzziness [default: near]

sampling type for fdr estimation of fuzziness change fdr. 
Either 'same', 'random' or 'near'.

#### -sn, --sample_names

names of samples seperated by ','. Example: 'sampleA,sampleB,sampleC'. If not provided, names are searched automaticly.

#### -cn, --control_names

names of controls seperated by ','. Example: 'controlA,controlB'. If not provided, names are searched automaticly.

### Interactive use

The script can be run interactively as well. Just provide the arguments in the script in the respective section on the top and source the script.

## Output:

Output are nucleosome positions files with new columns added for adjusted FDRs. All other columns stay the same. The files will be in a new directory which can be adjusted via the `-o, --output_dir` argument. (default: `/adjustedFDR/`)
