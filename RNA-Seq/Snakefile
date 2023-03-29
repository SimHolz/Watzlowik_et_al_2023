"""
Author: Simon Holzinger
Affiliation: AG Längst, University of Regensburg
Aim: Analysis of Plasmodium falciparum RNase-Seq Data
Date: 2021-06-09
Run: snakemake --use-conda Snakefile
Workflow:
    Preprocessing:
        - FastQC of raw reads
        - trimming of reads with Trimmomatic
        - FastQC of trimmed reads

"""

configfile: "config.yaml"

include: "1_code/0_smk/preprocessing.smk"
include: "1_code/0_smk/mapping.smk"

####################
# Main Rule
####################

rule all:
  input:
    "2_pipeline/00_qc/preprocessing/multiqc/preprocessing_multiqc_report.html", # preprocessing
    expand("2_pipeline/02_alignment/bigwigs/{sample}_{cases}.out.bw", sample = config["samples"], cases = ["Signal.Unique.str1","Signal.Unique.str2","Signal.UniqueMultiple.str1","Signal.UniqueMultiple.str2"]),
    "2_pipeline/00_qc/mapping/multiqc/mapping_multiqc_report.html",
