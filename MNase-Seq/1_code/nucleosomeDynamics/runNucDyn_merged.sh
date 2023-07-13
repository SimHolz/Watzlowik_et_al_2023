#!/bin/bash

# Bash Script for nucleosome dynamics package based on docker


# ADD ABSOLUTE PATH TO PROJECT HERE:
project_dir=<PATH TO PROJECT>



genome_annotation=/project_dir/0_data/reference/pf_nucDyn.gff
genome_sizes=/project_dir/0_data/reference/plasmoDB52/sizes.genome

##################
# merge bam files

mkdir -p ../../2_pipeline/05_NucleosomeDynamics/bam/
samtools merge ../../2_pipeline/05_NucleosomeDynamics/bam/iKO_24h.bam ../../2_pipeline/03_alignment_filtered/bam/iKO_24h_Rep1.s.bam ../../2_pipeline/03_alignment_filtered/bam/iKO_24h_Rep2.s.bam
samtools merge ../../2_pipeline/05_NucleosomeDynamics/bam/noKO_24h.bam ../../2_pipeline/03_alignment_filtered/bam/noKO_24h_Rep1.s.bam ../../2_pipeline/03_alignment_filtered/bam/noKO_24h_Rep2.s.bam

for SAMPLE in iKO_24h noKO_24h; do

	##################
	# sort bam files
	samtools sort -@ 8 -o ../../2_pipeline/05_NucleosomeDynamics/bam/$SAMPLE.s.bam ../../2_pipeline/05_NucleosomeDynamics/bam/$SAMPLE.bam

	# ##################
	# # remove unsorted bams
	rm ../../2_pipeline/05_NucleosomeDynamics/bam/$SAMPLE.bam

	##################
	# read bam
	mkdir -p ../../2_pipeline/05_NucleosomeDynamics/RData/
	docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited readBAM \
		--input /project_dir/2_pipeline/05_NucleosomeDynamics/bam/$SAMPLE.s.bam \
		--output /project_dir/2_pipeline/05_NucleosomeDynamics/RData/$SAMPLE.RData \
		--type paired

	##################
	# detect nucleosomes with nucleR
	mkdir -p ../../2_pipeline/05_NucleosomeDynamics/nucleR/
	docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited nucleR \
		--input /project_dir/2_pipeline/05_NucleosomeDynamics/RData/$SAMPLE.RData \
		--output /project_dir/2_pipeline/05_NucleosomeDynamics/nucleR/$SAMPLE.gff \
		--type paired --fragmentLen 175 --thresholdPercentage 50

	mkdir -p ../../2_pipeline/05_NucleosomeDynamics/statistics/nucleR
	docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited nucleR_stats \
		--input /project_dir/2_pipeline/05_NucleosomeDynamics/nucleR/$SAMPLE.gff \
		--genome $genome_annotation \
		--out_genes /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/nucleR/$SAMPLE''_genes.csv \
		--out_gw /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/nucleR/$SAMPLE''_genome.csv

	##################
	# detect Txstart
	mkdir -p ../../2_pipeline/05_NucleosomeDynamics/txstart/
	docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited txstart \
		--calls /project_dir/2_pipeline/05_NucleosomeDynamics/nucleR/$SAMPLE.gff \
		--genome $genome_annotation \
		--output /project_dir/2_pipeline/05_NucleosomeDynamics/txstart/$SAMPLE.gff

	mkdir -p ../../2_pipeline/05_NucleosomeDynamics/statistics/txstart/
	docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited txstart_stats \
		--input /project_dir/2_pipeline/05_NucleosomeDynamics/txstart/$SAMPLE.gff \
		--genome $genome_annotation \
		--out_genes /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/txstart/$SAMPLE''_genes.csv \
		--out_gw /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/txstart/$SAMPLE''_genome.png \
		--out_gw2 /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/txstart/$SAMPLE''_genome2.png
done

###################
# detect nucleosome dynamics
mkdir -p ../../2_pipeline/05_NucleosomeDynamics/NucDyn/
docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited nucDyn \
	--input2 /project_dir/2_pipeline/05_NucleosomeDynamics/RData/iKO_24h.RData \
	--input1 /project_dir/2_pipeline/05_NucleosomeDynamics/RData/noKO_24h.RData \
	--calls2 /project_dir/2_pipeline/05_NucleosomeDynamics/nucleR/iKO_24h.gff \
	--calls1 /project_dir/2_pipeline/05_NucleosomeDynamics/nucleR/noKO_24h.gff \
	--outputGff /project_dir/2_pipeline/05_NucleosomeDynamics/NucDyn/24h.gff \
	--outputBigWig /project_dir/2_pipeline/05_NucleosomeDynamics/NucDyn/24h.bw \
	--genome $genome_sizes

mkdir -p ../../2_pipeline/05_NucleosomeDynamics/statistics/NucDyn/
docker run -v $project_dir/:/project_dir/ nucleosomedynamics_docker_edited nucDyn_stats \
	--input /project_dir/2_pipeline/05_NucleosomeDynamics/NucDyn/24h.gff \
	--genome $genome_annotation \
	--out_genes /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/NucDyn/24h_genes.csv \
	--out_gw /project_dir/2_pipeline/05_NucleosomeDynamics/statistics/NucDyn/24h_genome.png
