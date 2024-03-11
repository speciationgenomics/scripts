#!/bin/bash

# Script to convert vcf to eigenstrat format for ADMIXTOOLS
# Written by Joana Meier
# It takes a single argument: the vcf file (can be gzipped) 
# optionally you can specify --renameScaff if you have scaffold names (not chr1, chr2...)
# The script will use a constant recombination rate of 2 cM/Mb assuming that you do not know the recombination rate. You may have to change that or change the jackknife block size in ADMIXTOOLS.

# Here, you can change the recombination rate which is currently set to 2 cM/Mb 
rec=2

# It requires vcftools and admixtools

# for some clusters, it is needed to load these modules:
# module load gcc/4.8.2 vcftools openblas/0.2.13_seq perl/5.18.4 admixtools

renameScaff="FALSE"

# If help is requested
if [[ $1 == "-h" ]]
then
 echo "Please provide the vcf file to parse, and optionally add --renameScaff if you have scaffolds instead of chromosomes"
 echo "Usage: convertVCFtoEigenstrat.sh <vcf file> --renameScaff --rec 2 (note, the second argument is optional)"
 exit 1

# If the second argument renameScaff is given, set it to True
elif [[ $2 == "--renameScaff" ]]
then
 renameScaff="TRUE"

# If no argument is given or the second one is not -removeChr, give information and quit
elif [ $# -ne 1 ]
then
 echo "Please provide the vcf file to parse, and optionally add --renameScaff if you have scaffolds instead of chromosomes"
 echo "Usage: ./convertVCFtoEigenstrat.sh <vcf file> --renameScaff (note, the second argument is optional)"
 exit 1
fi

# Set the first argument to the file name
file=$1
file=${file%.gz}
file=${file%.vcf}


# if the vcf file is gzipped:

if [ -s $file.vcf.gz ]
then

 # If renaming of scaffolds is requested, set all chromosome/scaffold names to 1
 if [ $renameScaff == "TRUE" ]
 then
  echo "setting scaffold names to 1 and positions to additive numbers"
  zcat $file".vcf.gz" | awk 'BEGIN {OFS = "\t";add=0;lastPos=0;scaff=""}{
    if($1!~/^#/){
       if($1!=scaff){add=lastPos;scaff=$1}
       $1="1"
       $2=$2+add
       lastPos=$2
     }
     print $0}' | gzip > $file.renamedScaff.vcf.gz

  # Get a .map and .ped file (remove multiallelic SNPs, monomorphic sites and indels)
  vcftools --gzvcf $file".renamedScaff.vcf.gz" \
         --plink --mac 1.0 --remove-indels --max-alleles 2 --out $file

 else
 # Get a .map and .ped file (remove multiallelic SNPs, monomorphic sites and indels)
 vcftools --gzvcf $file".vcf.gz" \
         --plink --mac 1.0 --remove-indels --max-alleles 2 --out $file
 fi

# if the file is not gzipped
else
 # If renaming of scaffolds is requested, set all chromosome/scaffold names to 1
 if [ $renameScaff == "TRUE" ]
 then
  echo "setting scaffold names to 1 and positions to additive numbers"
  awk 'BEGIN {OFS = "\t";add=0;lastPos=0;scaff=""}{
    if($1!~/^#/){
       if($1!=scaff){add=lastPos;scaff=$1}
       $1="1"
       $2=$2+add
       lastPos=$2
     }
     print $0}' $file.vcf | gzip > $file.renamedScaff.vcf.gz

  # Get a .map and .ped file (remove multiallelic SNPs, monomorphic sites and indels)
  vcftools --gzvcf $file".renamedScaff.vcf.gz" \
         --plink --mac 1.0 --remove-indels --max-alleles 2 --out $file
 else
 vcftools --vcf $file".vcf" --plink --mac 1.0 --remove-indels --max-alleles 2 --out $file
 fi
fi


# Change the .map file to match the requirements of ADMIXTOOLS by adding fake Morgan positions (assuming a recombination rate of 2 cM/Mbp)
awk -F"\t" -v rec=$rec 'BEGIN{scaff="";add=0}{
        split($2,newScaff,":")
        if(!match(newScaff[1],scaff)){
                scaff=newScaff[1]
                add=lastPos
        }
        pos=add+$4
	       count+=0.00000001*rec*(pos-lastPos)
        print newScaff[1]"\t"$2"\t"count"\t"pos
        lastPos=pos
}' ${file}.map  | sed 's/^chr//' > better.map
mv better.map ${file}.map

# Change the .ped file to match the ADMIXTOOLS requirements
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; 
 for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' ${file}.ped > tmp.ped
mv tmp.ped ${file}.ped


# create an inputfile for convertf
echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT.${file}
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT.${file}
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT.${file}
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT.${file}
echo "genotypeoutname: ${file}.eigenstratgeno" >> par.PED.EIGENSTRAT.${file}
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT.${file}
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT.${file}
echo "familynames:     NO" >> par.PED.EIGENSTRAT.${file}


# Use CONVERTF to parse PED to eigenstrat
convertf -p par.PED.EIGENSTRAT.${file}

# change the snp file for ADMIXTOOLS:
awk 'BEGIN{i=0}{i=i+1; print $1"\t"$2"\t"$3"\t"i"\t"$5"\t"$6}' $file.snp > $file.snp.tmp
mv $file.snp.tmp $file.snp
