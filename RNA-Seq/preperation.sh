#!/bin/bash

# Script for preperation of snakemake workflow

# Create reference directory
mkdir -p 0_data/reference/plasmoDB52/

# Load reference fasta and gff
wget -O 0_data/reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7_Genome.fasta https://plasmodb.org/common/downloads/release-52/Pfalciparum3D7/fasta/data/PlasmoDB-52_Pfalciparum3D7_Genome.fasta 
wget -O 0_data/reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7.gff https://plasmodb.org/common/downloads/release-52/Pfalciparum3D7/gff/data/PlasmoDB-52_Pfalciparum3D7.gff

# Create genome sizes file
samtools faidx 0_data/reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7_Genome.fasta
cut -f1,2 0_data/reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7_Genome.fasta.fai > 0_data/reference/plasmoDB52/sizes.genome