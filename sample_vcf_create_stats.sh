#!/bin/sh
# Mark Ravinet
# set variables
VCF=~/vcf/cichlid_full.vcf.gz
SUBSET_VCF=~/vcf/cichlid_subset.vcf.gz
OUT=~/vcftools/cichlid_subset

# randomly sample vcf and generate some statistics in order to set thresholds
bcftools view $VCF | vcfrandomsample -r 0.012 > $SUBSET_VCF

# allele freq
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --freq2 --out $OUT

# depth per ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --depth --out $OUT

# depth per site
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-mean-depth --out $OUT

# het per ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --het --out $OUT

# site quality
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-quality --out $OUT

# missing ind
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --missing-indv --out $OUT

# missing site
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --missing-site --out $OUT

# site quality
vcftools --gzvcf $SUBSET_VCF --min-alleles 2 --max-alleles 2 --site-quality --out $OUT
