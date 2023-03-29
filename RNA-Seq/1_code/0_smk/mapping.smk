"""
Author: Simon Holzinger
Affiliation: AG Längst, University of Regensburg
Aim: Mapping of Plasmodium RNA-Seq Data
Date: 2021-06-10
Workflow:
    - Create STAR Index
    - map reads with STAR
    - create bigwig tracks
    - move and index bams
    - QC
"""

####################
# Convert GFF to GTF
####################

rule convert_gff:
    input:
        gff=config["reference-annotation"]
    output:
        gtf="2_pipeline/00_annotation/annotation.gtf"
    conda:
        "../../5_envs/gffread.yaml"
    shell:
        "gffread {input.gff} -T -o {output.gtf}"

####################
# Create STAR Index
####################

rule index_STAR:
    input:
        fasta=config["reference-genome"],
        gtf="2_pipeline/00_annotation/annotation.gtf",
    output:
        "2_pipeline/00_index/STAR/Genome",
    params:
        basename="2_pipeline/00_index/STAR/",
        extra="--sjdbOverhang 56 "
            "--genomeSAindexNbases 11 "
    threads: config["parallel_threads"]
    benchmark:
        "2_pipeline/00_benchmarks/star-index/reference.star-index.benchmark.txt"
    conda:
        "../../5_envs/star.yaml"
    log:
        "2_pipeline/00_logs/star_index/reference.log"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.basename} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} {params.extra} > {log}"

####################
# Alignment STAR
####################

rule alignment_STAR:
    input:
        r1="2_pipeline/01_trimmed_reads/{sample}_1.fastq.gz",
        r2="2_pipeline/01_trimmed_reads/{sample}_2.fastq.gz",
        indexG="2_pipeline/00_index/STAR/Genome",
    output:
        temp("2_pipeline/02_alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam"),
        "2_pipeline/02_alignment/{sample}/{sample}_Log.final.out",
        temp("2_pipeline/02_alignment/{sample}/{sample}_Signal.Unique.str1.out.wig"),
        temp("2_pipeline/02_alignment/{sample}/{sample}_Signal.Unique.str2.out.wig"),
        temp("2_pipeline/02_alignment/{sample}/{sample}_Signal.UniqueMultiple.str1.out.wig"),
        temp("2_pipeline/02_alignment/{sample}/{sample}_Signal.UniqueMultiple.str2.out.wig"),
    params:
        index="2_pipeline/00_index/STAR/",
        basename="2_pipeline/02_alignment/{sample}/{sample}_",
        extra="--readFilesCommand gunzip -c "
            "--outSAMtype BAM SortedByCoordinate "
            "--outFilterType BySJout "
            "--outFilterMultimapNmax 20 "
            "--alignSJDBoverhangMin 1 "
            "--outFilterMismatchNmax 999 "
            "--outFilterMismatchNoverReadLmax 0.04 "
            "--alignIntronMin 5 "
            "--alignIntronMax 1000 "
            "--alignMatesGapMax 1000000 "
            "--outMultimapperOrder Random "
            "--outWigType wiggle ",
        seed=config["rndSeed"]
    threads: config["parallel_threads"]
    benchmark:
        "2_pipeline/00_benchmarks/alignment/{sample}.star.benchmark.txt"
    conda:
        "../../5_envs/star.yaml"
    log:
        "2_pipeline/00_logs/star/{sample}.log"
    shell:
        "STAR "
        "{params.extra} "
        "--runThreadN {threads} "
        "--genomeDir {params.index} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outFileNamePrefix {params.basename} "
        "--runRNGseed {params.seed} "
        "--outStd Log > {log}"

####################
# Move and Convert mapping Results
####################

rule wig_to_bigwig:
    input:
        "2_pipeline/02_alignment/{sample}/{sample}_{cases}.out.wig"
    output:
        "2_pipeline/02_alignment/bigwigs/{sample}_{cases}.out.bw"
    params:
        chrSize=config["chromSizes"],
    conda: "../../5_envs/wigtobigwig.yaml"
    shell:
        "wigToBigWig {input} {params.chrSize} {output}"

rule move_STAR:
    input:
        bam="2_pipeline/02_alignment/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "2_pipeline/02_alignment/{sample}.s.bam"
    shell:
        "cp {input} {output}"

rule index_bam:
    input:
        "2_pipeline/02_alignment/{sample}.s.bam"
    output:
        "2_pipeline/02_alignment/{sample}.s.bam.bai"
    shell:
        "samtools index {input}"

####################
# Mapping QC - Qualimap, Samtools
####################

rule qualimap:
    input:
        "2_pipeline/02_alignment/{sample}.s.bam"
    output:
        "2_pipeline/00_qc/mapping/qualimap/{sample}_qualimap/qualimapReport.html"
    params:
        annotationFile=config["reference-annotation"],
        outdir="2_pipeline/00_qc/mapping/qualimap/{sample}_qualimap/",
        extra="--java-mem-size=4G"
    conda:
        "../../5_envs/qualimap.yaml"
    log:
        "2_pipeline/00_logs/qualimap/{sample}.log"
    shell:
        "qualimap bamqc -bam {input} {params.extra} -gff {params.annotationFile} -outdir {params.outdir} -c > {log}"

rule samtools_flagstat:
    input:
        "2_pipeline/02_alignment/{sample}.s.bam"
    output:
        "2_pipeline/00_qc/mapping/samtools_flagstat/{sample}.txt"
    conda:
        "../../5_envs/samtools.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule samtools_idxstat:
    input:
        bam="2_pipeline/02_alignment/{sample}.s.bam",
        idx="2_pipeline/02_alignment/{sample}.s.bam.bai"
    output:
        "2_pipeline/00_qc/mapping/samtools_idxstats/{sample}_idxstats.txt"
    conda:
        "../../5_envs/samtools.yaml"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule samtools_stats:
    input:
        "2_pipeline/02_alignment/{sample}.s.bam"
    output:
        "2_pipeline/00_qc/mapping/samtools_stats/{sample}.txt"
    conda:
        "../../5_envs/samtools.yaml"
    shell:
        "samtools stats {input} > {output}"

####################
# Multiqc Report of alignment
####################
rule mapping_multiqc:
    input:
        expand("2_pipeline/00_qc/mapping/qualimap/{sample}_qualimap/qualimapReport.html", sample=config["samples"]),
        expand("2_pipeline/00_qc/mapping/samtools_flagstat/{sample}.txt", sample=config["samples"]),
        expand("2_pipeline/00_qc/mapping/samtools_idxstats/{sample}_idxstats.txt", sample=config["samples"]),
        expand("2_pipeline/00_qc/mapping/samtools_stats/{sample}.txt", sample=config["samples"])
    output:
        "2_pipeline/00_qc/mapping/multiqc/mapping_multiqc_report.html"
    params:
        star_dir="2_pipeline/02_alignment/",
        qualimap_dir="2_pipeline/00_qc/mapping/qualimap/",
        samtools_flagstat_dir="2_pipeline/00_qc/mapping/samtools_flagstat/",
        samtools_stats_dir="2_pipeline/00_qc/mapping/samtools_stats/",
        samtools_idxstats_dir="2_pipeline/00_qc/mapping/samtools_idxstats/",
        outdir="2_pipeline/00_qc/mapping/multiqc/",
        extra="-i mapping "
                "-f "
    conda:
        "../../5_envs/multiqc.yaml"
    shell:
        "multiqc {params.extra} -o {params.outdir} {params.star_dir} {params.qualimap_dir} {params.samtools_flagstat_dir} {params.samtools_stats_dir} {params.samtools_idxstats_dir}"
