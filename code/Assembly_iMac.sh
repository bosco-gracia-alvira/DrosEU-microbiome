#!/bin/bash
# This script assembles the microbial reads from each sample
# I am following the next Anvi'o post:
# https://merenlab.org/data/tara-oceans-mags/
# Bosco Gracia Alvira, 2023

### VARIABLES
CODE="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/code"
WORKDIR="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data"
RAW_READS="$WORKDIR"/raw_reads

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

### COMMANDS
if [[ ! -d "$WORKDIR"/DeNovo_Assembly ]]
then  
  mkdir -p "$WORKDIR"/DeNovo_Assembly
fi

cd "$WORKDIR"/DeNovo_Assembly || exit

anvi-run-workflow -w metagenomics \
                  -c "$CODE"/metagenomics-config.json \
                  --save-workflow-graph




if [[ ! -d "$WORKDIR" ]]
then  
  mkdir -p \
    "$WORKDIR"/00_LOGS \
    "$WORKDIR"/02_ASSEMBLY \
    "$WORKDIR"/03_CONTIGS \
    "$WORKDIR"/04_MAPPING \
    "$WORKDIR"/05_ANVIO_PROFILE
fi

if [[ ! -d "$TEMP" ]]
then  
  mkdir "$TEMP"
fi

if [[ ! -f "$TEMP"/Total_1.fq.gz ]]
then
  cat "$RAW_READS"/extracted/*_1.fastq.gz > "$TEMP"/Total_1.fq.gz
  cat "$RAW_READS"/extracted/*_2.fastq.gz > "$TEMP"/Total_2.fq.gz
fi

#Now we are ready for assembly. I will use SPAdes, that is among the most popular assemblers and gave the best results to Xiaomeng.
#https://cab.spbu.ru/files/release3.15.2/manual.html
echo "Assembling the metagenome with metaSPAdes"
metaspades.py \
  -1 "$TEMP"/Total_1.fq.gz \
  -2 "$TEMP"/Total_2.fq.gz \
  -o "$WORKDIR"/02_ASSEMBLY/ \
  -k 21,33,55,77 \
  -t 16

echo "Assembly done"

#This command renames the contigs and removes those smaller than 1 kb.
echo "Renaming contigs and removing those smaller than 1 kb"
anvi-script-reformat-fasta \
  "$WORKDIR"/02_ASSEMBLY/contigs.fasta \
  -o "$WORKDIR"/03_CONTIGS/contigs.fa \
  --min-len 1000 \
  --simplify-names \
  --report "$WORKDIR"/00_LOGS/name_conversions.txt

#With the following commands I map the pool-seq reads individually against the consensus contigs file. anvi-init-bam sorts and indexes bam files.
echo "Mapping the reads from each poolseq against the consensus contigs files"
bowtie2-build "$WORKDIR"/03_CONTIGS/contigs.fa "$WORKDIR"/04_MAPPING/contigs
for i in ${UNIQGEN};
do
  bowtie2 \
    --threads 8 \
    -x "$WORKDIR"/04_MAPPING/contigs \
    -1 "$READSDIR"/${i}/Pooled_${i}_Clean_1.fq.gz \
    -2 "$READSDIR"/${i}/Pooled_${i}_Clean_2.fq.gz \
    -S "$WORKDIR"/04_MAPPING/${i}.sam \
    -p 16 \
    2>"$WORKDIR"/00_LOGS/bowtie2_${i}.txt;
  samtools view -F 4 -bS "$WORKDIR"/04_MAPPING/${i}.sam > "$WORKDIR"/04_MAPPING/${i}-RAW.bam;
  anvi-init-bam "$WORKDIR"/04_MAPPING/${i}-RAW.bam -o "$WORKDIR"/04_MAPPING/${i}.bam;
  rm "$WORKDIR"/04_MAPPING/${i}.sam "$WORKDIR"/04_MAPPING/${i}-RAW.bam;
done
rm "$WORKDIR"/04_MAPPING/*.bt2
rm -r "$TEMP"
echo "Voilà, your metagenome is ready for binning"
