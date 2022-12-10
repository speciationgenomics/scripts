#!/bin/bash

# Script by Henry North

for (( i = 1; i < 101; i++ ))
do

  echo "running simulation ${i} of 100"

  ### 1: Run the simulation
  slim neutral_stable.slim > stable_${i}.out # run the simulation
  tail -n+15 stable_${i}.out | bgzip -c > stable_${i}.vcf.gz # zip and tidy up the vcf output
  rm stable_${i}.out # remove a temporary file (raw SLIM output)

  ### 2: Calculate summary statistics

  ## Calculate Fst from the VCF

  vcftools --gzvcf stable_${i}.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out stable_fst_${i} # produces a log file
  grep "weighted" stable_fst_${i}.log | awk '{print $7}' > fst_stable_${i}.txt # Extract the fst value from the vcftools log file
  rm stable_fst_${i}.weir.fst # remove a temporary file (genome-wide FST; we don't need this)
  rm stable_fst_${i}.log # remove a temporary file (the vcftools log file)


  ## Calculate nucleotide diversity from the VCF

  vcftools --gzvcf stable_${i}.vcf.gz --site-pi --out stable_pi_${i}
  awk 'NR > 1 {sum+=$3} END {print sum / (NR - 1)}' stable_pi_${i}.sites.pi > stable_pi_mean_${i}.txt  # calculate the mean
  rm stable_pi_${i}.log # remove intermediate logfile
  rm stable_pi_${i}.sites.pi # Remove the intermediate file

  ## Calculate Tajima's D from the VCF
  vcftools --gzvcf stable_${i}.vcf.gz --out stable_${i} --TajimaD 10000
  rm stable_${i}.log # remove the intermediate logfile
  awk 'NR==2' stable_${i}.Tajima.D | cut -f 4 > stable_tajD_${i}.txt # extract out the value of D
  rm stable_${i}.Tajima.D # remove the intermediate file

  rm stable_${i}.vcf.gz # remove a temporary file (vcf; we don't need this any more)

done


### Collate Fst data
cat fst_stable_*.txt >> all_fst_stable.txt # Put all 100 fst values into one text output
rm fst_stable_*.txt # remove the individual fst text files

# As above for pi and D
cat stable_pi_mean_*.txt >> all_pi_stable.txt
rm stable_pi_mean_*.txt

cat stable_tajD_*.txt >> all_tajD_stable.txt
rm stable_tajD_*.txt
