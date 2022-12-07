for (( i = 1; i < 101; i++ )); do
  slim neutral_expansion.slim > neutral_${i}.out # this took a solid 5 minutes for one run

  tail -n+15 neutral_${i}.out | bgzip -c > neutral_${i}.vcf.gz
  rm neutral_${i}.out

  vcftools --gzvcf neutral_${i}.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out neutral_${i}
  rm neutral_${i}.vcf.gz

  grep "weighted" neutral_1.log | awk '{print $7}' > fst_neutral_${i}.txt
  rm neutral_${i}.weir.fst
done
