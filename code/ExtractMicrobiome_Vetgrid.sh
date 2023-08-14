#!/bin/bash
# This script maps each poolseq against the reference D. melanogaster and D. simulans genomes and retains only the unmapped reads.
# Part of this script was taken from:
# https://github.com/capoony/DrosEU_pipeline/tree/master
# Bosco Gracia Alvira, 2023

### VARIABLES
RAW_READS="/Volumes/Temp/DrosEU"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data"
REFERENCE="$WORKDIR/references/reference"

### COMMANDS

# We check if there are raw reads and if we have a reference folder
if [[ ! -d "$RAW_READS"/extracted ]]
then
    mkdir "$RAW_READS"/extracted
fi

if [[ ! -d "$WORKDIR"/references ]]
then
    mkdir -p \
        "$WORKDIR"/references
fi

if [[ ! -d "$WORKDIR"/raw_reads ]]
then
    mkdir -p \
        "$WORKDIR"/raw_reads
fi

# We create a reference genome by concatenating D. melanogaster (v6.52) and D. simulans (v2.02) genomes from flybase.org, as well as other eukaryotic genomes that could potentially be contaminating our samples based on previous experiences.
# These eukaryotes are: H. sapiens, M. musculus, A. thaliana, S. cerevisiae and C. lupus familiaris. We download the latest reference genome from RefSeq.
if [[ ! -f "$WORKDIR"/references/reference.fa ]]
then
    # D. simulans
    wget -O "$WORKDIR"/references/dsim202.fa.gz \
        http://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.02_FB2020_03/fasta/dsim-all-chromosome-r2.02.fasta.gz
    
    # D. melanogaster
    wget -O "$WORKDIR"/references/dmel652.fa.gz \
        http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.52_FB2023_03/fasta/dmel-all-chromosome-r6.52.fasta.gz
    
    # H. sapiens
    wget -O "$WORKDIR"/references/Hsapiens.fa.gz \
        ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

    # M. musculus
    wget -O "$WORKDIR"/references/Mmusculus.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/reference/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz

    # A. thaliana
    wget -O "$WORKDIR"/references/Athaliana.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_thaliana/reference/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

    # S. cerevisiae
    wget -O "$WORKDIR"/references/Scerevisiae.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

    # C. lupus
    wget -O "$WORKDIR"/references/Clupus.fa.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/representative/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz

    # We concatenate the 7 genomes into a single file and index it.
    gzcat "$WORKDIR"/references/dsim202.fa.gz \
          "$WORKDIR"/references/dmel652.fa.gz \
          "$WORKDIR"/references/Hsapiens.fa.gz \
          "$WORKDIR"/references/Mmusculus.fa.gz \
          "$WORKDIR"/references/Athaliana.fa.gz \
          "$WORKDIR"/references/Scerevisiae.fa.gz\
          "$WORKDIR"/references/Clupus.fa.gz > "$WORKDIR"/references/reference.fa

    # We index the reference
    bwa index -p "$WORKDIR"/references/reference "$WORKDIR"/references/reference.fa

    # We don't want the genomes wasting space anymore
    rm "$WORKDIR"/references/*.fa.gz

fi


# Reads retrieved from ENA are already trimmed and don't have adapters :) If you don't believe me run FastQC
# We align them to the reference file

for i in $(cut -f1 "$WORKDIR"/ENA_* | grep -v "run_accession")
do
    if [[ -f "$RAW_READS"/"${i}"_1.fastq.gz ]]
    then
        bwa mem \
            -M \
            -t 24 \
            -R "@RG\tID:${i}\tSM:${i}" \
            "$REFERENCE" \
            "$RAW_READS"/"${i}"_1.fastq.gz \
            "$RAW_READS"/"${i}"_2.fastq.gz |\
        samtools view \
            -Sbh -F 0x100 -@ 24 -f 0xc - |\
        samtools collate -Ou -@ 24 - |\
        samtools fixmate -u -@ 24 - - |\
        samtools view -u -@ 24 -f 0x1 - |\
        samtools fastq -@ 24 -N -0 /dev/null -s /dev/null \
            -1 "$WORKDIR"/raw_reads/"${i}"_1.fastq.gz \
            -2 "$WORKDIR"/raw_reads/"${i}"_2.fastq.gz \
            -
    elif [[ ! -f "$RAW_READS"/"${i}"_1.fastq.gz ]]
    then
        bwa mem \
            -M \
            -t 24 \
            -R "@RG\tID:${i}\tSM:${i}" \
            "$REFERENCE" \
            "$RAW_READS"/"${i}".fastq.gz |\
        samtools view \
            -Sbh -F 0x100 -f 0x4 - |\
        samtools collate -Ou -@ 24 - |\
        samtools fastq -@ 24 -0 /dev/null - > "$WORKDIR"/raw_reads/"${i}".fastq.gz
    else
        echo "Sample ${i} could not be processed" >> "$WORKDIR"/raw_reads/log.txt
    fi
done

# We check if the same accessions are found in the "raw reads" and in the "extracted reads"
for i in $(cut -f1 "$WORKDIR"/ENA_* | grep -v "run_accession")
do
    if [[ -n "$(find "$RAW_READS" -name "${i}*.fastq.gz" -maxdepth 1 -print -quit)" && -n "$(find ""$WORKDIR"/raw_reads" -name "${i}*.fastq.gz" -print -quit)" ]]
    then
        :
    else
        echo "Sample ${i} is missing" >> "$WORKDIR"/raw_reads/log2.txt
    fi
done

# Using the ENA information we rename the files to make them match the library name
for i in $(cat "$WORKDIR"/ENA_* | cut -f1,2 | tr "\t" ",")
do  
    rename -s $(echo "${i}" | cut -f1 -d ",") $(echo "${i}" | cut -f2 -d ",") "$WORKDIR"/raw_reads/"$(echo "${i}" | cut -f1 -d ",")"_?.fastq.gz
done
