#!/bin/bash
# This script downloads D. melanogaster poolseq data from BioProjects PRJNA388788 and PRJNA308584. These data were published in the following papers, respectively: 
# https://doi.org/10.1093/molbev/msaa120
# https://doi.org/10.7554/eLife.67577
# Bosco Gracia Alvira, 2023

### VARIABLES
RAW_READS="/Volumes/Temp/DrosEU"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data"
READSDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Ancestral_microbiome/data/R1_mapping/00_READS"
METADATA=""$WORKDIR"/metadata.csv"
REFERENCE="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads/dsimM252v1.2+microbiome.fa"

### COMMANDS
if [[ ! -d "$RAW_READS" ]]
then  
  mkdir -p \
    "$RAW_READS"
fi

esearch -db sra -query 'PRJNA388788[bioproject]' | efetch -format runinfo > "$WORKDIR"/PRJNA388788.csv

esearch -db sra -query 'PRJNA308584[bioproject]' | efetch -format runinfo > "$WORKDIR"/PRJNA308584.csv



prefetch SRR8494439 | fasterq-dump SRR8494439

