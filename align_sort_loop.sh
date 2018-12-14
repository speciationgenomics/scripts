#!/bin/sh
# Mark Ravinet

# a simple script to map lots of individuals
REF=~/reference/P_nyererei_v2.fasta

INDS=($(for i in ~/raw/*R1.100k.fastq.gz; do echo $(basename ${i%.R*}); done))

for IND in ${INDS[@]};
do
	# declare variables
	FORWARD=~/raw/${IND}.R1.100k.fastq.gz
	REVERSE=~/raw/${IND}.R2.100k.fastq.gz
	OUTPUT=~/align/${IND}_sort.bam

	# then align and sort
	echo "Aligning $IND with bwa"
	bwa mem -M -t 2 $REF $FORWARD \
	$REVERSE | samtools view -b | \
	samtools sort -T ${IND} > $OUTPUT

done
