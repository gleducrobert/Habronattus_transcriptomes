from __future__ import division
import io
import os
from collections import defaultdict
#from itertools import izip_longest
from itertools import izip_longest as iz
import itertools
from os import listdir
from os.path import isfile, join
import sys

#argument 1 = directory to go through
#argument 2 = output directory name

directory = sys.argv[1]
outputDirectory = sys.argv[2]

d = os.path.dirname(outputDirectory)
if not os.path.exists(outputDirectory):
	os.makedirs(outputDirectory)

onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f))]

for f in onlyfiles:
	if f.endswith(".fa") or f.endswith(".fas") or f.endswith(".fasta"):
	
		print("reading file " + f )

		lines = io.open(directory + "/" + f, 'r', newline="") 

		outfile = {}

		for line in lines:
			line = line.rstrip('\n')
			line = line.rstrip('\r')
			
			if line.startswith(">"):			
				line = line[1:]
				info = line.split('|')
				species = info[0]
				if species not in outfile:
					outfile[species] = ">" + species + '\n'
	
			else:
				sequence = line.rstrip("\n")
				outfile[species] = outfile[species] + sequence + "\n"
		

		fileName = f 
		outputFile = io.open(outputDirectory + "/" + fileName, 'w') 
		
		for species in outfile:
			outputFile.write(outfile[species])
