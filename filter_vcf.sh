#!/bin/sh
# Mark Ravinet

VCF_IN=~/vcf/cichlid_full.vcf.gz
VCF_OUT=~/vcf/cichlid_full_filtered.vcf.gz

# set filters
MAF=0.1
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50

# filter vcf for analysis
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > \
$VCF_OUT
