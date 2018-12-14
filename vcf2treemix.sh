#!/bin/bash

# Script to convert vcf file to Treemix format
# by Joana Meier
# requires vcftools, plink and plink2treemix.py from Treemix package to be installed
# takes as input the vcf file and the clust file

if [ $# -ne 2 ]
 then
 echo "Please provide the following arguments: <vcf file> <clust file>"
 echo "The .clust file contains three columns: samplename\tsamplename\tgroup"
 exit 1
fi

clust=$2
file=${1%.gz}
file=${file%.vcf}

# Use VCFtools to make a map and a ped file, using only bi-allelic SNPs with mac 2 (also creates a log file)
if [ -s $file.vcf.gz ]
then

 # Get a .map and .ped file
 vcftools --gzvcf $file".vcf.gz" \
         --plink --mac 2 --remove-indels --max-alleles 2 \
         --out $file

else
 file=${file%.vcf}
 vcftools --vcf $file".vcf" \
         --plink --mac 2 --remove-indels --max-alleles 2  \
         --out $file

fi

# Adjust the map file to allow for non-human chromosome names (else the new plink version struggles)
awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
        print $0
}' ${file}.map > better.map
mv better.map ${file}.map

# convert it to a stratified frq file, also creates .bed, .bim, .fam, .log, .nosex
plink --file $file --make-bed --out $file --allow-no-sex --allow-extra-chr 0
plink --bfile $file --freq --missing --within $2 --out $file --allow-no-sex --allow-extra-chr 0

# zip it
gzip $file".frq.strat"

# create input file for Treemix
plink2treemix.py $file".frq.strat.gz" $file".treemix.frq.gz"

# unzip allele frequency information
gunzip $file".treemix.frq.gz"
gunzip $file".frq.strat.gz"

# make a file with the positions
awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' $file".map" > $file".positions"
paste $file".positions" $file".treemix.frq" > $file".frequencies"

awk '{printf $0
	for(i = 4; i <= NF; i++){
		split($i,values,",")
		if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2])
		else freq=0
		printf freq"\t"
	}
	printf "\n"}' $file".frequencies" > $file".frequencies2"
mv $file".frequencies2" $file".frequencies"

awk 'BEGIN{scaffold="";pos=0;newpos=0}
	{if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' \
	$file".frequencies" > $file".frequencies.newpos"

gzip $file".treemix.frq"
