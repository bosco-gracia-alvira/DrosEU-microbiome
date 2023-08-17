#!/bin/bash
# This script assembles the microbial reads from each sample
# I am following the next Anvi'o post:
# https://merenlab.org/data/tara-oceans-mags/
# Bosco Gracia Alvira, 2023

### VARIABLES
CODE="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/code"
WORKDIR="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data/DeNovo_Assembly"
RAW_READS="/Volumes/Data/Dropbox (PopGen)/Bosco/PhD_Dropbox/DrosEU-microbiome/data/raw_reads"

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

### COMMANDS
if [[ ! -d "$WORKDIR" ]]
then  
  mkdir -p "$WORKDIR"
fi

cd "$WORKDIR" || exit

# Generate samples.txt

# First I remove any previous tmp files
rm "$WORKDIR"/*.tmp

# Anvi'o requires a samples.txt file with four columns: sample name, group (samples in the same group will be co-assembled), and paths to forward and reverse reads
echo -e "sample\tgroup\tr1\tr2" > "$WORKDIR"/header.tmp
for i in $(basename "$RAW_READS"/*_1.fastq.gz)
do
  echo "${i%_1.fastq.gz}" >> "$WORKDIR"/sample.tmp
  echo "DrosEU" >> "$WORKDIR"/group.tmp
  echo "../raw_reads/$i" >> "$WORKDIR"/r1.tmp
  echo "../raw_reads/${i%_1.fastq.gz}_2.fastq.gz" >> "$WORKDIR"/r2.tmp
done

paste \
  "$WORKDIR"/sample.tmp \
  "$WORKDIR"/group.tmp \
  "$WORKDIR"/r1.tmp \
  "$WORKDIR"/r2.tmp > "$WORKDIR"/body.tmp

cat "$WORKDIR"/header.tmp "$WORKDIR"/body.tmp > "$WORKDIR"/samples.txt

rm "$WORKDIR"/*.tmp

anvi-run-workflow -w metagenomics \
                  -c "$CODE"/metagenomics-config.json

# In case it is needed
anvi-run-workflow -w metagenomics \
                  -c "$CODE"/metagenomics-config.json \
                  --additional-params \
                    --unlock \
                    --keep-going \
                    --rerun-incomplete

echo "Voilà, your metagenome is ready for binning"
