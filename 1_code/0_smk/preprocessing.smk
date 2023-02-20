"""
Author: Simon Holzinger
Affiliation: AG Längst, University of Regensburg
Aim: Preprocessing of Plasmodium RNA-Seq Data
Date: 2021-06-10
Workflow:
    - FastQC of raw reads
    - trimming of reads with Trimmomatic
    - FastQC of trimmed reads
    - Create Report with Multiqc
"""


####################
# Raw FastQC
####################

rule raw_fastqc:
    input:
        lambda wildcards: f"{config['samples'][wildcards.sample]}_{wildcards.num}_001.fastq.gz"
    output:
        html="2_pipeline/00_qc/preprocessing/raw_fastqc/{sample}_{num}_fastqc.html",
        zip="2_pipeline/00_qc/preprocessing/raw_fastqc/{sample}_{num}_fastqc.zip"
    wrapper:
        "0.35.0/bio/fastqc"

####################
# Trimming of Reads
####################

rule trim_reads:
    input:
        r1= lambda wildcards: f"{config['samples'][wildcards.sample]}_R1_001.fastq.gz",
        r2= lambda wildcards: f"{config['samples'][wildcards.sample]}_R2_001.fastq.gz"
    output:
        r1="2_pipeline/01_trimmed_reads/{sample}_1.fastq.gz",
        r2="2_pipeline/01_trimmed_reads/{sample}_2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="2_pipeline/01_trimmed_reads/{sample}_1.unpaired.fastq.gz",
        r2_unpaired="2_pipeline/01_trimmed_reads/{sample}_2.unpaired.fastq.gz",
        log="2_pipeline/00_logs/trimmomatic/{sample}.log"
    log:
        "2_pipeline/00_logs/trimmomatic/{sample}.log"
    params:
        adapt_trim="ILLUMINACLIP:" + config["adapters"] + ":2:30:10",
        qual_trim="MAXINFO:30:0.2",
        len_trim="MINLEN:35"
    threads: 12
    benchmark:
        "2_pipeline/00_benchmarks/trimmomatic/{sample}.trimmomatic.benchmark.txt"
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} {params.adapt_trim} {params.qual_trim} {params.len_trim} 2> {log}"


####################
# FastQC of trimmed reads
####################
rule trimmed_fastqc:
    input:
        "2_pipeline/01_trimmed_reads/{sample}_{n}.fastq.gz"
    output:
        html="2_pipeline/00_qc/preprocessing/trimmed_fastqc/{sample}_{n}_fastqc.html",
        zip="2_pipeline/00_qc/preprocessing/trimmed_fastqc/{sample}_{n}_fastqc.zip"
    params:
        outdir="2_pipeline/00_qc/trimmed_fastqc/"
    log:
        "2_pipeline/00_logs/trimmed_fastqc/{sample}_{n}.log"
    wrapper:
        "0.35.0/bio/fastqc"


####################
# Multiqc Report of preprocessing
####################
rule preprocessing_multiqc:
    input:
        expand("2_pipeline/00_qc/preprocessing/raw_fastqc/{sample}_{num}_fastqc.html", sample=config["samples"], num = ["R1", "R2"]),
        expand("2_pipeline/00_qc/preprocessing/trimmed_fastqc/{sample}_{n}_fastqc.html", sample=config["samples"], n = ["1","2"]),
        expand("2_pipeline/00_logs/trimmomatic/{sample}.log", sample=config["samples"])
    output:
        "2_pipeline/00_qc/preprocessing/multiqc/preprocessing_multiqc_report.html"
    params:
        raw_fastqc_dir="2_pipeline/00_qc/preprocessing/raw_fastqc/",
        trimmomatic_dir="2_pipeline/00_logs/trimmomatic/",
        trimmed_fastqc_dir="2_pipeline/00_qc/preprocessing/trimmed_fastqc/",
        outdir="2_pipeline/00_qc/preprocessing/multiqc/",
        extra="-i preprocessing "
                " -f "
                " -d "
                " -dd 1 "
    conda:
        "../../5_envs/multiqc.yaml"
    shell:
        "multiqc {params.extra} -o {params.outdir} {params.raw_fastqc_dir} {params.trimmomatic_dir} {params.trimmed_fastqc_dir}"
