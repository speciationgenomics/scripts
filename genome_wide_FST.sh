##!/bin/sh
# Mark Ravinet


VCF1=~/vcf/cichlid_full_filtered.vcf.gz
VCF2=~/vcf/cichlid_full_filtered_rename.vcf.gz

# rename individuals
bcftools query -l $VCF1 > samples
#
# remove paths
while read ind
do
  basename $ind
done < samples > samples2

# create sample names
bcftools reheader -s samples2 -o $VCF2 $VCF1

# index
bcftools index $VCF2

# make the population files
grep "PunPundMak" samples2 > ppmak
grep "PunPundPyt" samples2 > pppyt
grep "PunNyerPyt" samples2 > pnpyt
grep "PunNyerMak" samples2 > pnmak

# calculate genoem wide FastQC
vcftools --gzvcf ${VCF2} \
--weir-fst-pop ppmak --weir-fst-pop pppyt --out ../vcftools/ppmak_ppyt
vcftools --gzvcf ${VCF2} \
--weir-fst-pop ppmak --weir-fst-pop pnpyt --out ../vcftools/ppmak_pnpyt
vcftools --gzvcf ${VCF2} \
--weir-fst-pop ppmak --weir-fst-pop pnmak --out ../vcftools/ppmak_pnmak
vcftools --gzvcf ${VCF2} \
--weir-fst-pop pppyt --weir-fst-pop pnpyt --out ../vcftools/pppyt_pnpyt
vcftools --gzvcf ${VCF2} \
--weir-fst-pop pnpyt --weir-fst-pop pnmak --out ../vcftools/pnpyt_pnmak
