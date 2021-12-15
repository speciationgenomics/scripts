#!/usr/bin/env python2

# read in gzip file

import gzip
import sys, argparse
import glob, re
import numpy as np

def getOptionValue(option): # needs sys
  if option in sys.argv:
    optionPos = sys.argv.index(option)
    optionValue = sys.argv[optionPos + 1]
    return optionValue
  else:
    print >> sys.stderr, "\nWarning, option", option, "not_specified.\n"

if "-i" in sys.argv:
  fileName = getOptionValue("-i")
else:
  print "\nplease specify input file name using -i <file_name> \n"
  sys.exit()
  
if "-o" in sys.argv:
  prefix = getOptionValue("-o")
else:
  print "\nplease specify output prefix using -o \n"
  sys.exit()

# prefix = "test"
# fileName = "./test.ld.gz"

f = gzip.open(fileName, 'r')

f.readline()

counter = 0
chr_distance_ld = {}

for line in f:

	line = line.strip().split()
	distance = abs(int(line[1])-int(line[4]))
	chr = line[0]		
	if chr in chr_distance_ld:
		if distance in chr_distance_ld[chr]:
			chr_distance_ld[chr][distance].append(float(line[6]))
		else:
			chr_distance_ld[chr][distance] = [float(line[6])]
	else:
		chr_distance_ld[chr] = {}
		chr_distance_ld[chr][distance] = [float(line[6])]
		#print (chr_distance_ld)		
	counter += 1		
	if (counter % 100000000) == 0:
		print ("Read %d entires." % counter)
		sys.stdout.flush()
		
		
print("Calculating average distances.")	
sys.stdout.flush()

decay_output = open(prefix + ".ld_decay", "w")
#write header
decay_output.write("chr\tdistance\tavg_R2\tstd\n")	
#decay_output.write("distance\tR2\n")
decay_output_bins = open(prefix + ".ld_decay_bins", "w")
#write header
decay_output_bins.write("chr\tdistance\tavg_R2\tstd\n")	
#decay_output.write("distance\tR2\n")	
	
#distance_avg_ld = {}
for chr in chr_distance_ld:
	distance_bin = {}
	for distance in sorted(chr_distance_ld[chr]):
		lds = np.array(chr_distance_ld[chr][distance])
		mean = np.mean(lds)
		std = np.std(lds)
		#sum(chr_distance_ld[chr][distance])/len(chr_distance_ld[chr][distance]
		decay_output.write("%s\t%d\t%g\t%g\n" % (chr, distance, mean, std))
		#round to closest 1000, and add to new structure No, round to 500 (all from 0-1000).
		bin_thousand = int(round(distance-500, -3)) + 500
		if bin_thousand in distance_bin:
			distance_bin[bin_thousand] = np.concatenate([distance_bin[bin_thousand], lds])
		else:
			distance_bin[bin_thousand] = lds
	for bin in sorted(distance_bin):
		mean = np.mean(distance_bin[bin])
		std = np.std(distance_bin[bin])
		decay_output_bins.write("%s\t%d\t%g\t%g\n" % (chr, bin, mean, std))
		
decay_output_bins.close()
decay_output.close()

f.close()
