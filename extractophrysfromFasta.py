#This script is to extract the transcript of a specific species from each contig file, and rename
#the fasta header so that it is the name of the contig (it can be the file name, then I can remove the extra parts
#of the name using textwrangler.

from os import listdir
from os.path import isfile, join
import os
import itertools
import io
import sys

directory = sys.argv[1]

onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f))]
ophrys_file = open(sys.argv[2], "w")
outfile = {}

for f in onlyfiles:

	if f.endswith(".fa") or f.endswith(".fas") or f.endswith(".fasta"):
		print("About to read file " + f)
		 
		lines = io.open(directory + "/" + f, 'r', newline="") 
		contigName = f
		species_sequence = {}
		
		for line in lines:
			if line.startswith(">"):
				line = line[1:] 
				words = line.rstrip('\n').split("|")
				species = words[0] 
				length = words[1]
				if species == "ophrys":
					outfile[contigName] = ">" + str(contigName) + "|" + str(length) + "|\r"
					print(contigName)
			else:
				sequence = line
				species_sequence[species] = sequence
		if contigName in outfile:
				outfile[contigName] = outfile[contigName] + species_sequence["ophrys"] + "\r"

for contigName in outfile:
	ophrys_file.write(outfile[contigName])
