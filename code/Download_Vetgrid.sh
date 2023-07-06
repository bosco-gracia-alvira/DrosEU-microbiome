#!/bin/bash
# This script downloads D. melanogaster poolseq data from BioProjects PRJNA388788 and PRJNA308584. These data were published in the following papers, respectively: 
# https://doi.org/10.1093/molbev/msaa120
# https://doi.org/10.7554/eLife.67577
# Bosco Gracia Alvira, 2023

### VARIABLES
RAW_READS="/Volumes/Temp/DrosEU"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data"

### COMMANDS

# We check if the folder where the raw reads will be downloaded exists
if [[ ! -d "$RAW_READS" ]]
then
  mkdir -p \
    "$RAW_READS"
fi

# We check if we have the ENA files needed to download the reads. 
# These files have to be downloaded manually from the ENA webpage.
# I select the fields: "run_accession", "experimental_alias", "fastq_ftp" and "library_name".
if [[ ! -f "$WORKDIR"/ENA_* ]]
then
  echo "You need at least one ENA file if you want to download reads."
  exit
fi

# We create a txt file with the ftp paths. The forth column contains the ftp paths.
cat "$WORKDIR"/ENA_* | cut -f4 | tr ";" "\n" | grep -v "fastq_ftp" > "$WORKDIR"/ftp.txt

# We download all the files into "RAW_READS"
wget -P "$RAW_READS" -i "$WORKDIR"/ftp.txt

# Using the ENA information we rename the files to make them match the library name
for i in $(cat "$WORKDIR"/ENA_* | cut -f1,2 | tr "\t" ",")
do  rename "s/$(echo "${i}" | cut -f1 -d ",")\./$(echo "${i}" | cut -f2 -d ",")\./" "$RAW_READS"/*
done

