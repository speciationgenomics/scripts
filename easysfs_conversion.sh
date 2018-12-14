#!/bin/sh
# Mark Ravinet

# convert to SFS
# set variables
VCF=~/vcf/chr20_filtered.vcf.gz

# make pops
bcftools query -l $VCF | grep "Mak" | awk '{split($0,a,"."); print $1,a[2]}' > pop_file

# use this line to estimate projections
easySFS.py -i $VCF -p pop_file -a -f --preview

# this line will calculate the necessary SFS files
easySFS.py -i $VCF -p pop_file -a -f --proj 8,8
