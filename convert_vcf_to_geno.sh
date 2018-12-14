#!/bin//sh
# Mark Ravinet

# convert vcf to geno
VCF=~/cichlids.vcf.gz
OUTFILE=~/${CHR}.geno


# for sparrows
# initialise file
echo -ne "scaffold\tpos\t" > $OUTFILE
# echo sample names and also deal with running line
bcftools query -l $VCF | tr "\n" "\t" >> $OUTFILE
echo "" >> $OUTFILE
# output IUPAC genotypes and remove non calls
bcftools query -e 'STRLEN(REF)>1' -f '%CHROM\t%POS[\t%IUPACGT]\n' $VCF | sed -e 's/\.\/\./N/g' >> $OUTFILE
gzip $OUTFILE
