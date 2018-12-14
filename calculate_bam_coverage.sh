#!/bin/sh
# Mark Ravinet

module load samtools

for BAM in *_sort.bam
do
  # calculate mean and standard deviation of bam coverage
  STATS=$(samtools depth -a $BAM | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}')
  echo -e "${BAM}\t${STATS}"
done > bam_coverage.txt
