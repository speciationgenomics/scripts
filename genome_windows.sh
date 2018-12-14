#!/bin/sh
# Mark Ravinet

### genome wide estimates

INPUT=~/genome_scan/chr20.geno.gz
OUTPUT1=~/genome_scan/Pundamilia.ABBABABA.w20k.s20k.csv
OUTPUT2=~/genome_scan/Pundamilia.withKivu.dxy.pi.csv.gz

# calculate ABBABABA for NyerMak introgression into nyererei-like at Python (for files with SNPs only, thus m=10)
ABBABABAwindows.py -w 20000 -m 10 -s 20000 -g $INPUT -o $OUTPUT1 \
   -f phased -T 5 --minData 0.5 --writeFailedWindows\
   -P1 PundPyt 11725.PunPundPyt,11727.PunPundPyt,11728.PunPundPyt,11729.PunPundPyt \
   -P2 NyerPyt 11719.PunNyerPyt,11986.PunNyerPyt,11992.PunNyerPyt,11546.PunNyerPyt \
   -P3 NyerMak 11591.PunNyerMak,11593.PunNyerMak,11595.PunNyerMak,11598.PunNyerMak \
   -O kivu 64253 \

# calculate Dxy and pi including the outgroup to compare divergence between the species to divergence between Victoria and Kivu
popgenWindows.py -w 20000 -m 10000 -s 20000 -g $INPUT -o $OUTPUT2 \
   -f phased -T 5 \
   -p PundPyt 11725.PunPundPyt,11727.PunPundPyt,11728.PunPundPyt,11729.PunPundPyt \
   -p NyerPyt 11719.PunNyerPyt,11986.PunNyerPyt,11992.PunNyerPyt,11546.PunNyerPyt \
   -p NyerMak 11591.PunNyerMak,11593.PunNyerMak,11595.PunNyerMak,11598.PunNyerMak \
   -p PundMak 13069.PunPundMak,10558.PunPundMak,10560.PunPundMak,11297.PunPundMak \
   -p kivu 64253
