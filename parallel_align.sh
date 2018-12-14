#!/bin/sh
# Mark Ravinet

# align parallel
REF=~/reference/P_nyererei_v2.fasta

echo "Aligning $1 with bwa"
bwa mem -M -t 4 $REF ./data/${1}.R1.fastq.gz \
./data/${1}.R2.fastq.gz | samtools view -b  > ./align/${1}.bam
