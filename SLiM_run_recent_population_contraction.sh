#!/bin/bash

# Script by Henry North

for (( i = 1; i < 101; i++ ))
do

  echo "running simulation ${i} of 100"

  ### 1: Run the simulation
  slim neutral_recentContraction.slim > recentContraction_${i}.out # run the simulation
  tail -n+15 recentContraction_${i}.out | bgzip -c > recentContraction_${i}.vcf.gz # zip and tidy up the vcf output
  rm recentContraction_${i}.out # remove a temporary file (raw SLIM output)

  ### 2: Calculate summary statistics

  ## Calculate Fst from the VCF

  vcftools --gzvcf recentContraction_${i}.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out recentContraction_fst_${i} # produces a
log file
  grep "weighted" recentContraction_fst_${i}.log | awk '{print $7}' > fst_recentContraction_${i}.txt # Extract the fst value from the vcftools
log file
  rm recentContraction_fst_${i}.weir.fst # remove a temporary file (genome-wide FST; we don't need this)
  rm recentContraction_fst_${i}.log # remove a temporary file (the vcftools log file)


  ## Calculate nucleotide diversity from the VCF

  vcftools --gzvcf recentContraction_${i}.vcf.gz --site-pi --out recentContraction_pi_${i}
  awk 'NR > 1 {sum+=$3} END {print sum / (NR - 1)}' recentContraction_pi_${i}.sites.pi > recentContraction_pi_mean_${i}.txt  # calculate the
mean
  rm recentContraction_pi_${i}.log # remove intermediate logfile
  rm recentContraction_pi_${i}.sites.pi # Remove the intermediate file

  ## Calculate Tajima's D from the VCF
  vcftools --gzvcf recentContraction_${i}.vcf.gz --out recentContraction_${i} --TajimaD 10000
  rm recentContraction_${i}.log # remove the intermediate logfile
  awk 'NR==2' recentContraction_${i}.Tajima.D | cut -f 4 > recentContraction_tajD_${i}.txt # extract out the value of D
  rm recentContraction_${i}.Tajima.D # remove the intermediate file

  rm recentContraction_${i}.vcf.gz # remove a temporary file (vcf; we don't need this any more)

done


### Collate Fst data
cat fst_recentContraction_*.txt >> all_fst_recentContraction.txt # Put all 100 fst values into one text output
rm fst_recentContraction_*.txt # remove the individual fst text files

# As above for pi and D
cat recentContraction_pi_mean_*.txt >> all_pi_recentContraction.txt
rm recentContraction_pi_mean_*.txt

cat recentContraction_tajD_*.txt >> all_tajD_recentContraction.txt
rm recentContraction_tajD_*.txt
