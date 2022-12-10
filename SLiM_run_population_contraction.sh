#!/bin/bash

# Script by Henry North

for (( i = 1; i < 101; i++ ))
do

  echo "running simulation ${i} of 100"

  ### 1: Run the simulation
  slim neutral_contraction.slim > contraction_${i}.out # run the simulation
  tail -n+15 contraction_${i}.out | bgzip -c > contraction_${i}.vcf.gz # zip and tidy up the vcf output
  rm contraction_${i}.out # remove a temporary file (raw SLIM output)

  ### 2: Calculate summary statistics

  ## Calculate Fst from the VCF

  vcftools --gzvcf contraction_${i}.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out contraction_fst_${i} # produces a log file
  grep "weighted" contraction_fst_${i}.log | awk '{print $7}' > fst_contraction_${i}.txt # Extract the fst value from the vcftools log file
  rm contraction_fst_${i}.weir.fst # remove a temporary file (genome-wide FST; we don't need this)
  rm contraction_fst_${i}.log # remove a temporary file (the vcftools log file)


  ## Calculate nucleotide diversity from the VCF

  vcftools --gzvcf contraction_${i}.vcf.gz --site-pi --out contraction_pi_${i}
  awk 'NR > 1 {sum+=$3} END {print sum / (NR - 1)}' contraction_pi_${i}.sites.pi > contraction_pi_mean_${i}.txt  # calculate the mean
  rm contraction_pi_${i}.log # remove intermediate logfile
  rm contraction_pi_${i}.sites.pi # Remove the intermediate file

  ## Calculate Tajima's D from the VCF
  vcftools --gzvcf contraction_${i}.vcf.gz --out contraction_${i} --TajimaD 10000
  rm contraction_${i}.log # remove the intermediate logfile
  awk 'NR==2' contraction_${i}.Tajima.D | cut -f 4 > contraction_tajD_${i}.txt # extract out the value of D
  rm contraction_${i}.Tajima.D # remove the intermediate file

  rm contraction_${i}.vcf.gz # remove a temporary file (vcf; we don't need this any more)

done


### Collate Fst data
cat fst_contraction_*.txt >> all_fst_contraction.txt # Put all 100 fst values into one text output
rm fst_contraction_*.txt # remove the individual fst text files

# As above for pi and D
cat contraction_pi_mean_*.txt >> all_pi_contraction.txt
rm contraction_pi_mean_*.txt

cat contraction_tajD_*.txt >> all_tajD_contraction.txt
rm contraction_tajD_*.txt
