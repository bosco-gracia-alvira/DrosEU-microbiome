#!/bin/bash
# This script bins the DrosEU metagenomic assembly obtained in the previous script and imports them to the Anvio db.
# I am following the these posts:
# https://merenlab.org/data/tara-oceans-mags/
# https://thegeneschool.github.io/metagenomics/anvio/ (only for generating the bins collection.txt)
# Bosco Gracia Alvira, 2023

### VARIABLES
CODE="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/code"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data/DeNovo_Assembly"
RAW_READS="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data/raw_reads"

### COMMANDS
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

# Security check...
if [[ ! -d "$WORKDIR"/05_ANVIO_PROFILE ]]
then 
  echo "You don't have the folder 05_ANVIO_PROFILE, maybe you should run Assembly_Vetgrid.sh."
  exit
else
  if [[ -d "$WORKDIR"/06_MERGED ]]
  then  
    mkdir "$WORKDIR"/07_BINNING
  else
    echo "You don't have the folder 06_MERGED, which means that you haven't finished running the Assembly_Vetgrid_wf.sh script. Maybe you should..."
    exit
  fi
fi

# In principle I use MetaBat2 as binner

if [[ ! -d "$WORKDIR"/07_BINNING/MetaBat2 ]]
then  
  mkdir "$WORKDIR"/07_BINNING/MetaBat2
fi

jgi_summarize_bam_contig_depths \
  --outputDepth "$WORKDIR"/07_BINNING/MetaBat2/Metabat2depth.txt \
  "$WORKDIR"/04_MAPPING/DrosEU/*.bam

metabat2 \
  --seed 999 \
  -i "$WORKDIR"/02_FASTA/DrosEU/DrosEU-contigs-prefix-formatted-only.fa \
  -a "$WORKDIR"/07_BINNING/MetaBat2/Metabat2depth.txt \
  -m 2500 \
  -o "$WORKDIR"/07_BINNING/MetaBat2/METABAT_ \
  -t 16

rename 's/_\./_/' "$WORKDIR"/07_BINNING/MetaBat2/METABAT*.fa

cd "$WORKDIR"/07_BINNING/metabat2

# This perl command that I have definitely not created myself creates a two-columns file that links each contig of the metagenome to the bin in which it has been allocated
grep ">" METABAT*.fa | perl -p -i -e "s/METABAT(.+).fa:>(.+)/\$2\tMETABAT\$1/g" > collection.txt

# We input the collection into Anvi'o
anvi-import-collection \
    "$WORKDIR"/07_BINNING/MetaBat2/collection.txt \
    --contigs-mode \
    -p "$WORKDIR"/06_MERGED/DrosEU/PROFILE.db \
    -c "$WORKDIR"/03_CONTIGS/DrosEU-contigs.db \
    -C METABAT2

# We summarise the collection
anvi-summarize \
    -p "$WORKDIR"/06_MERGED/DrosEU/PROFILE.db \
    -c "$WORKDIR"/03_CONTIGS/DrosEU-contigs.db \
    -C METABAT2 \
    -o "$WORKDIR"/07_BINNING/METABAT_SUMMARY

anvi-refine \
    -p "$WORKDIR"/06_MERGED/DrosEU/PROFILE.db \
    -c "$WORKDIR"/03_CONTIGS/DrosEU-contigs.db \
    -C METABAT2 \
    -b METABAT_119
