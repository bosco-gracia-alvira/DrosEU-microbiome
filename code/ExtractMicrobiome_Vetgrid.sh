#!/bin/bash
# This script maps each poolseq against the reference D. melanogaster and D. simulans genomes and retains only the unmapped reads.
# Part of this script was taken from:
# https://github.com/capoony/DrosEU_pipeline/tree/master
# Bosco Gracia Alvira, 2023

### VARIABLES
RAW_READS="/Volumes/Temp/DrosEU"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data"
REFERENCE="$WORKDIR"/references/reference.fa.gz
ADAPTERS="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/db/Adapters"
### COMMANDS

# We check if there are raw reads and if we have a reference folder
if [[ ! -f "$RAW_READS"/*.fastq.gz ]]
then
    echo "You don't have any poolseq from which you can extract the microbiome."
    exit
fi

if [[ ! -d "$WORKDIR"/references ]]
then
    mkdir -p \
        "$WORKDIR"/references
fi

# We create a reference genome by concatenating D. melanogaster (v6.52) and D. simulans (v2.02) genomes from flybase.org
wget -O "$WORKDIR"/references/dsim202.fa.gz http://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.02_FB2020_03/fasta/dsim-all-chromosome-r2.02.fasta.gz
wget -O "$WORKDIR"/references/dmel652.fa.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.52_FB2023_03/fasta/dmel-all-chromosome-r6.52.fasta.gz

zcat "$WORKDIR"/references/dmel652.fa.gz "$WORKDIR"/references/dsim202.fa.gz > "$WORKDIR"/references/reference.fa.gz


# We use bbduk, from bbtools, to trim the reads. 
# The first command removes the illumina adapters. The adapters.fa file can be downloaded from the BBMap documents.
for i in $(basename "$RAW_READS"/*_1.fastq.gz | cut -f1 -d "_")
do  
    bbduk.sh \
        -Xmx24g \
        in1="$RAW_READS"/${i}_1.fastq.gz \
        in2="$RAW_READS"/${i}_2.fastq.gz \
        out1="$RAW_READS"/${i}.RmAdp_1.fq.gz \
        out2="$RAW_READS"/${i}.RmAdp_2.fq.gz \
        ref="$ADAPTERS"/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe
done

#We asses the quality of the reads after removing the adapters.
for i in $(cat $RAW/Sample.list)
do
fastqc \
    $RM/fastq_wo_adapt/${i}.RmAdp_1.fq.gz \
    $RM/fastq_wo_adapt/${i}.RmAdp_2.fq.gz \
    -o $RM/fastqc_wo_adapt
done
#We use bbduk again, now for trimming any PhiX174 sequence, that is used in Illumina sequencing and remove the low quality reads (Q-score < 20). Again, phix174_ill.ref.fa was downloaded from BBMap. More info about BBDuk:
#https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/
for i in $(cat $RAW/Sample.list)
do
bbduk.sh \
    in1=$RM/fastq_wo_adapt/${i}.RmAdp_1.fq.gz \
    in2=$RM/fastq_wo_adapt/${i}.RmAdp_2.fq.gz \
    out1=$RM/fastq_clean/${i}.clean_1.fq.gz \
    out2=$RM/fastq_clean/${i}.clean_2.fq.gz \
    ref=~/PhD/db/Adapters/phix174_ill.ref.fa \
    k=31 hdist=1 qtrim=rl trimq=20
done

#We asses the quality of the reads after quality trimming. If you look at the html files you will realise that the quality of the reverse reads increases a lot.
for i in $(cat $RAW/Sample.list)
do
fastqc \
    $RM/fastq_clean/${i}.clean_1.fq.gz \
    $RM/fastq_clean/${i}.clean_2.fq.gz \
    -o $RM/fastqc_clean
done



# Trim reads 
cutadapt \
-q 18 \
--minimum-length 75 \
-o trimmed-read1.fq.gz \
-p trimmed-read2.fq.gz \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
-B CAAGCAGAAGACGGCATACGAGAT \
-O 15 \
-n 3 \
read1.fq.gz \
read2.fq.gz




# align to reference
bwa mem \
-M \
-t 24 \
"$REFERENCE" \
trimmed-read1.fastq.gz \
trimmed-read2.fastq.gz \
| samtools view \
-Sbh -q 20 -F 0x100 - > library.bam

# create reads groups 