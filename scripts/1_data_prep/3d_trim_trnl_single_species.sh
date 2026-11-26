# Note: run from base directory of repository
#       for example: ./scripts/1_data_prep/3d_trim_trnl_single_species.sh 

#!/bin/bash

# Define inpout and output directory
INDIR="data/sequences/trnl/single_species/raw_reads/"
OUTDIR="data/sequences/trnl/single_species/trimmed/"

# Create output directory
mkdir $OUTDIR

# Trim sequences in each file
for file in "${INDIR}"*R1_001.fastq.gz
do
    #  Get sample name from file path
    SAMPLE=$(echo ${file} | sed "s/_R1_001\.fastq\.gz//")
    # Status update
    echo ${SAMPLE##*/}_R1_001.fastq.gz ${SAMPLE##*/}_R2_001.fastq.gz
    # Trim random Illumina adapter (5 bases) and primer
    # Linking primer and reverse complement of paired primer
    # 1. Specify forward linked adapter including random sequence (5 bp)
    # 2. Specify reverse linked adapter including random sequenc (5 bp)
    # 3. Minimum length = 10 maximum length = 100 -> discard untrimmed
    # 3 cont. Minimum match of adapte/primer combo = 17 bases
    # 4. Specify paird out files
    # 5. Specify paired in files
    cutadapt -a ^N{5}GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGG \
        -A ^N{5}CCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC \
        -m 10 -M 100 --discard-untrimmed -O 17 \
        -o ${SAMPLE}_R1.fastq -p ${SAMPLE}_R2.fastq \
        ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz
done
