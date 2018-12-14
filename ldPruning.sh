#!/bin/bash

# script that prunes SNPs with too high linkage 
# requires vcftools

#@author: Joana Meier
#@date: December, 2016

# Check if help is requested
if [[ $1 = "-h" ]]; then
        echo -e "\nUsage: ldPruning.sh <vcffile[.gz]> [optional: <LD threshold (R^2), default 0.1> <outputformat vcf/plink> <if 0/1 recoding requested, write 01>]"
  exit 0
fi

# Read in the file name (remove potential file endings given)
file=$1
file=${file%.gz}
file=${file%.vcf}

# Set default values
thresh=0.1
format="vcf"

# If the vcf file is zipped:
gz1=""
gz2=""
if [[ $1 == *"vcf.gz" ]]; then
  gz1="gz"
  gz2=".gz"
fi


# Check if the provided file exists and contains data
if [[ ! -s $file".vcf" ]] && [[ ! -s $file".vcf.gz" ]]
then
        echo -e "\n"$1" not found or empty. Please provide a vcf file with data"
        echo -e "\nUsage: ldPruning.sh <vcffile[.gz]> [optional: <LD threshold (R^2), default 0.1> <outputformat vcf/plink> <if 0/1 recoding requested, write TRUE>]"
        exit 1
fi


# Check if more than one argument is provided which ones
if (( "$#" > 1 ))    
then
        case $2 in
	  01)
		r="01"
		;;
	  vcf)
		format="vcf"
		;;
	  plink)
		format="plink"
		;;
	  *)
		thresh=$2
		;;
	esac
fi

if (( "$#" > 2 ))
then
        case $3 in
          01)
                r="01"
                ;;
          vcf)
                format="vcf"
                ;;
          plink)
                format="plink"
                ;;
          *)
                thresh=$3
                ;;
	esac
fi 

if (( "$#" > 3 ))
then
        case $4 in
          01)
                r="01"
                ;;
          vcf)
                format="vcf"
                ;;
          plink)
                format="plink"
                ;;
          *)
                thresh=$4
                ;;
	esac
fi


if (( "$#" > 4 ))
then
        echo "Error: Too many arguments provided"
        echo -e "\nUsage: ldPruning.sh <vcffile> [optional: <LD threshold (R^2), default 0.1> <outputformat vcf/plink> <if 0/1 recoding requested, write 01>]"
        exit 1
fi

# Let the user know, that it is working on it
echo "working..."

# For running on the Euler cluster, load the required modules
#module load plink/1.07
#module load openblas/0.2.13_par
#module load zlib/1.2.8
#module load vcftools


# Check which SNPs are in too high linkage and output a list of SNPs to be pruned out
vcftools --${gz1}vcf ${file}.vcf${gz2} --plink --out ${file} 2> tmp

# Set all chromosomes to 1 to trick plink into accepting non-human chr names
awk '{split($2,chr,":"); $1=1; $2="1:"chr[2]; print $0}' ${file}.map > ${file}.map.tmp
mv ${file}.map.tmp ${file}.map
plink --file $file --indep-pairwise 50 10 $thresh --out $file --noweb --silent

sed -i 's/:/\t/g' ${file}.prune.in


# Output pruned file either in vcf or plink format (if requested, 01 recoded) 
if (( format == "vcf" ))
then
	vcftools --${gz1}vcf ${file}.vcf${gz2} --out $file.pruning --positions $file.prune.in --stdout --recode > $file.LDpruned.vcf
else 
	vcftools --${gz1}vcf ${file}.vcf${gz2} --positions $file.prune.in --out $file.LDpruned --plink
	if (( r == "01" ))
	then
		plink --file $file.LDpruned  --recode$r --out $file.LDpruned$r 
	fi
fi

# Clean up unnecessary info file of vcftools which does not want to be silent at all 
rm tmp

# Output info about number of pruned SNPs to the console
echo "finished, new file "$file.LDpruned$r" filtered for LD in 50 kb windows, shifting by 10 kb with LD threshold "$thresh
echo `grep "After filtering" $file.pruning.log`
