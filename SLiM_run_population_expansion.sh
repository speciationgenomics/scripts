#!/bin/bash 

# Script by Henry North

for (( i = 1; i < 101; i++ ))
do

  echo "running simulation ${i} of 100"

  ### 1: Run the simulation
  slim neutral_expansion.slim > expansion_${i}.out # run the simulation
  tail -n+15 expansion_${i}.out | bgzip -c > expansion_${i}.vcf.gz # zip and tidy up the vcf output
  rm expansion_${i}.out # remove a temporary file (raw SLIM output)

  ### 2: Calculate summary statistics

  ## Calculate Fst from the VCF

  vcftools --gzvcf expansion_${i}.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out expansion_fst_${i} # produces a log file
  grep "weighted" expansion_fst_${i}.log | awk '{print $7}' > fst_expansion_${i}.txt # Extract the fst value from the vcftools log file
  rm expansion_fst_${i}.weir.fst # remove a temporary file (genome-wide FST; we don't need this)
  rm expansion_fst_${i}.log # remove a temporary file (the vcftools log file)


  ## Calculate nucleotide diversity from the VCF

  vcftools --gzvcf expansion_${i}.vcf.gz --site-pi --out expansion_pi_${i}
  awk 'NR > 1 {sum+=$3} END {print sum / (NR - 1)}' expansion_pi_${i}.sites.pi > expansion_pi_mean_${i}.txt  # calculate the mean
  rm expansion_pi_${i}.log # remove intermediate logfile
  rm expansion_pi_${i}.sites.pi # Remove the intermediate file

  ## Calculate Tajima's D from the VCF
  vcftools --gzvcf expansion_${i}.vcf.gz --out expansion_${i} --TajimaD 10000
  rm expansion_${i}.log # remove the intermediate logfile
  awk 'NR==2' expansion_${i}.Tajima.D | cut -f 4 > expansion_tajD_${i}.txt # extract out the value of D
  rm expansion_${i}.Tajima.D # remove the intermediate file

  rm expansion_${i}.vcf.gz # remove a temporary file (vcf; we don't need this any more)

done


### Collate Fst data
cat fst_expansion_*.txt >> all_fst_expansion.txt # Put all 100 fst values into one text output
rm fst_expansion_*.txt # remove the individual fst text files

# As above for pi and D
cat expansion_pi_mean_*.txt >> all_pi_expansion.txt
rm expansion_pi_mean_*.txt

cat expansion_tajD_*.txt >> all_tajD_expansion.txt
rm expansion_tajD_*.txt
