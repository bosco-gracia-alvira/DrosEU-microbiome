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
for i in "$WORKDIR"/ENA_*.tsv
do
  if [[ ! -e "$i" ]]
  then
    echo "You need at least one ENA file if you want to download reads.
          The nomenclature I use is ENA_PRJNAXXXXXX.tsv"
    exit
  fi
done

# We create a txt file with the ftp paths. The forth column contains the ftp paths.
cat "$WORKDIR"/ENA_* | cut -f4 | tr ";" "\n" | grep -v "fastq_ftp" > "$WORKDIR"/ftp.txt

# We download all the files into "RAW_READS"
wget -P "$RAW_READS" -i "$WORKDIR"/ftp.txt

# We want to make sure that we have downloaded both pairs for all the files. 
# If it is not the case, we re-run wget on the lines that are missing.
for i in $(cut -f1 "$WORKDIR"/ENA_* | grep -v "run_accession")
do
    if [ -f "$RAW_READS"/"${i}"_1.fastq.gz ] && [ ! -f "$RAW_READS"/"${i}"_2.fastq.gz ]
    then
      echo "${i}_2.fastq.gz"
      wget -P "$RAW_READS" $(grep "${i}_2.fastq.gz" "$WORKDIR"/ftp.txt)
    elif [ -f "$RAW_READS"/"${i}"_2.fastq.gz ] && [ ! -f "$RAW_READS"/"${i}"_1.fastq.gz ]
    then 
      echo "${i}_1.fastq.gz"
      wget -P "$RAW_READS" $(grep "${i}_1.fastq.gz" "$WORKDIR"/ftp.txt)
    elif [ ! -f "$RAW_READS"/"${i}"_2.fastq.gz ] && [ ! -f "$RAW_READS"/"${i}"_1.fastq.gz ] && [ ! -f "$RAW_READS"/"${i}".fastq.gz ]
    then
      wget -P "$RAW_READS" $(grep "${i}" "$WORKDIR"/ftp.txt)
    fi
done

# For some reason, in paired-end SRR accessions a file with low quality reads is also downloaded. We don't want it.
# for i in $(cut -f1 "$WORKDIR"/ENA_* | grep -v "run_accession")
# do
#     if [ -f "$RAW_READS"/"${i}"_1.fastq.gz ] && [ -f "$RAW_READS"/"${i}"_2.fastq.gz ] && [ -f "$RAW_READS"/"${i}".fastq.gz ]
#     then
#       #rm "$RAW_READS"/"${i}".fastq.gz
#       echo "File "${i}".fastq.gz was lost forever. It weighted $(wc -c "$RAW_READS"/"${i}".fastq.gz | cut -f2 -d " ")"
#     fi
# done

