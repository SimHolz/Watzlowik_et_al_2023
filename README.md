# Watzlowik et al paper

This repository contains supplementary code for the Plasmodium Snf2L Paper. The subdirectories have a more specific readme file for installation, usage, etc.

## RNA-Seq

The `RNA-Seq` folder contains code used in the analysis of RNA-Seq data in the paper.
It has the following content:

- Snakemake pipeline for RNA-Seq preprocessing and mapping.
- Test data and preperation script to test the pipeline
- R Script for analyzing RNA-Seq data with DESeq2 (Sample metadata file required and not provided) and basic GO enrichment analysis.
- Non comprehensive R Script with code for mimicking plots shown in the publication.
- environment files for conda to use within the snakemake pipeline (`--use-conda`).

## MNase-Seq

The `MNase-Seq` folder contains code used in the analysis of RNA-Seq data in the paper.
It has the following content:

- Snakemake pipeline for MNase-Seq preprocessing, mapping, filtering and nucleosome position calling
- Test data and preperation script to test the pipeline
- environment files for conda to use within the snakemake pipeline (`--use-conda`)
- shell script to run the nucleosome dynamics package from a docker instance
- zipped docker file for nucleosome dynamics - test data removed for size reasons
- AdjustFDR Script to adjust DANPOS FDRs
- `find_dynamic_nucleosomes.R` a script to extract nucleosomes which show big differences in occupancy, fuzziness or shifts from DANPOS results.
- Non comprehensive R Script with code for mimicking plots shown in the publication.
