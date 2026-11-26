# Note: run from base directory of repository
#       for example: ./scripts/1_data_prep/3a_trim_its2_edna.sh 

#!/bin/bash

# Define inpout and output directory
INDIR="data/sequences/its2/edna/raw_reads/"
OUTDIR="data/sequences/its2/edna/trimmed/"

# Create output directory
mkdir $OUTDIR

# Trim sequences in each file
for file in "${INDIR}"*R1_001.fastq.gz
do
    #  Get sample name from file path
    SAMPLE=$(echo ${file##*/} | sed "s/_R1_001\.fastq\.gz//")
    # Status update
    echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz
    # Trim random Illumina adapter (5 bases) and primer
    # As the read lenth is shorter than the expected amplicon we're not using linked primers
    # 1. Specify forward linked adapter including random sequence (5 bp)
    # 2. Specify reverse linked adapter including random sequenc (5 bp)
    # 3. Minimum length = 15 (no maximum length as the expected lenght > 250 BP exeeds read lenght)
    # 3 cont. Minimum match of adapte/primer combo = 17 bases
    # 4. Specify paird out files
    # 5. Specify paired in files
    cutadapt -g ^N{5}ATGCGATACTTGGTGTGAAT \
        -G ^N{5}TCCTCCGCTTATTGATATGC \
        -m 15 --discard-untrimmed -O 20 \
        -o ${OUTDIR}${SAMPLE}_R1.fastq -p ${OUTDIR}${SAMPLE}_R2.fastq \
        ${INDIR}${SAMPLE}_R1_001.fastq.gz ${INDIR}${SAMPLE}_R2_001.fastq.gz
done

