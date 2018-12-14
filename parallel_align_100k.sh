#!/bin/sh
# Mark Ravinet

# align parallel
REF=~/reference/P_nyererei_v2.fasta

echo "Aligning $1 with bwa"
bwa mem -M -t 4 $REF ./data_100kreads/${1}.R1.100k.fastq.gz \
./data_100kreads/${1}.R2.100k.fastq.gz | samtools view -b | \
samtools sort -T ${1} > ./align_100kreads/${1}_sort.bam
