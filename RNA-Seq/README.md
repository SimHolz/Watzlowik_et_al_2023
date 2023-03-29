# RNA-Seq - Watzlowik et al

This directory contains code to reproduce the RNA-Seq analysis for the Plasmodium Snf2L Paper.
It contains the pipeline used to map the raw sequencing reads and a small test dataset to test the pipeline. Addionally there is example code for the analysis and plotting of the results.

## Pipeline

### Requirements

This code has been tested on macOS 13.4.1 but should run on any UNIX system.

It requires the following software:

  - conda (tested with v23.1.0)
  - samtools (tested with v1.12)
  - snakemake (tested with v7.26.0)
    - snakemake will install:
      - gffread v0.12.1
      - multiqc v1.11
      - qualimap v2.2.2d
      - samtools v1.12
      - star v2.7.9a
      - trimmomatic v0.39
      - ucsc-wigtobigwig v277

### Installation

Download the repository and install [conda](https://conda.io/projects/conda/en/stable/user-guide/install/index.html). Conda will handle installation of additional packages.

### Usage

Run the following commands to create a new environment with the required software.

```shell
conda create -n watzlowik_2023_RNA -c bioconda snakemake=7.26.0 samtools=1.12
```

Then activate the environment and make sure you are located in the RNA-Seq directory.
Now run the `preperation.sh` script to download reference genome and annotation.

```shell
conda activate watzlowik_2023_RNA
cd <PATH TO RNA-Seq DIR>
./preperation.sh
```

Finally adjust the `config.yaml` to your needs or just start the pipeline with the test-data provided.

```shell
snakemake --use-conda
```

To deactivate the conda environment run:

```shell
conda deactivate
```

### Expected Results

The snakemake pipeline creates an addional directory called `2_pipeline/` which contains logs for the individual steps (`00_logs`), quality control which is conveniently summarised as a multiqc report (`00_qc/preprocessing/multiqc/` and `00_qc/mapping/multiqc/`), intermediate trimmed fastq files (`01_trimmed_reads`) and resulting alignments (`02_alignment`) in `.bam` and `.bigwig` format

## Additional Code

Additional R-Markdown-code is provided in `0_code/1_r/` to show how analysis was done. However execution is limited with the provided test data. Although the scripts have been edited and documented to the best of our knowledge and belief the scripts have to be seen just as guidance and might not work immediately as is. Paths may require adjustments but the specific commands are the ones used in our analysis.

`DE_Analysis.Rmd` contains code to replicate the differential expression analysis done in this publication.

`plots.Rmd` contains code to recreate the plots in the manuscript.

### Requirements

Analysis was performed using R v4.2.2.
Additional libraries used are:

  - tidyverse v2.0.0
  - tximport v1.26.1
  - GenomicFeatures v1.50.4
  - DESeq2 v1.38.3
  - Rsubread v2.12.3
  - org.Pf.plasmo.db v3.12.0
  - topGO v2.50.0
  - DEGreport v1.34.0
  - ggpubr v0.6.0
