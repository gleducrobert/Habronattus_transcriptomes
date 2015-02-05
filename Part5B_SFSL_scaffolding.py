from __future__ import division
import io
import os
from collections import defaultdict
import itertools

numContigPerSpecies = {}

outfile = {}
#Species as the key and the sequence (and fused sequences) as the values. This will eventually be the new contig file.
SpeciesSE = {}
#dictionary of lists, where the key is the species and the value is a list of starts and ends
SpeciesN = {}
#dictionary of lists, where the key is the species and the value is a list of strings of Ns to add to sequences
species_sequence = {}
#dictionary of lists of species sequences, where the key is the species and the value is a list of sequence fragments
species_fused_sequence = {}
#dictionary of items of species fused sequences, where the key is the species and the value is the scaffolded transcript

species_fused_sequence_B = defaultdict(list)

outputFile = open("Pcontig_comp148941_c0_seq1.fa_B", "w")

fasta = io.open('Pcontig_comp148941_c0_seq1.fa', 'r', newline="") 

sequence = ''

for line in fasta:
	line = line.rstrip('\n')
	line = line.rstrip('\r')
	if line.startswith(">"):
		
##PART 1: SORTING THE SPECIES AND HEADER INFO
			
		line = line[1:]
		info = line.split('|')
		species = info[0]
		
		if species not in numContigPerSpecies:
			numContigPerSpecies[species] = 1
		else:
			numContigPerSpecies[species] += 1
		
		contig = info[1]
		start = (int(info[2]))
		end = (int(info[3]))
			
		#start by adding the species, and starting a blank list for each, to keep track of starts and end positions for each species
		if species not in SpeciesSE:
			list_of_positions = SpeciesSE.get(species, [])
		
		if species not in species_sequence:
			sequence_list = species_sequence.get(species,[]) 
		
		list_of_positions.append(start)
		list_of_positions.append(end)	
		list_of_positions.sort()		
		SpeciesSE[species] = list_of_positions
			
		if species not in outfile:
			outfile[species] = ">" + species + "|" + contig + "|" + str(start) + "|" 
		elif species in outfile:
			outfile[species] = outfile[species] + str(end) + "|"
#This is adding the header info to the new file content.

	else: 
	
#PART 2: SORTING THE SEQUENCE LINES
		sequence = line.rstrip("\n")
		sequence_list.append(sequence)
		species_sequence[species] = sequence_list

#now let's verify that the length(sequence_list) = numContigperSpecies[species]

for species in species_sequence:
	number_seq = len(species_sequence[species])
	print("number of sequences for " + species + " is " + str(number_seq))
	if number_seq != numContigPerSpecies[species]:
		print(number_seq)
		print(species + " has bad sequences")
	elif number_seq == numContigPerSpecies[species]:
		print(species + " has good sequences")

#GOOD! The sequence list matches the number of contigs per species.

#PART 3: SORTING THE NUMBER OF MISSING SEQUENCES AND ADDING THESE TO THE SEQUENCES
for species in SpeciesSE:
	list_of_positions = SpeciesSE[species]
	length = len(list_of_positions)
#	if length > 2:
	list_of_Ns = SpeciesN.get(species, [])
	if length == 2:
		print(species, length)
		SpeciesN[species] = list_of_Ns
		list_of_Ns.append("")
		print(list_of_Ns)
	if length == 4:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		list_of_Ns.append("N"*number_missing)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 6:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 8:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		number_missing_3 = list_of_positions[6] - list_of_positions[5]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		list_of_Ns.append("N"*number_missing_3)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 10:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		number_missing_3 = list_of_positions[6] - list_of_positions[5]
		number_missing_4 = list_of_positions[8] - list_of_positions[7]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		list_of_Ns.append("N"*number_missing_3)
		list_of_Ns.append("N"*number_missing_4)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 12:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		number_missing_3 = list_of_positions[6] - list_of_positions[5]
		number_missing_4 = list_of_positions[8] - list_of_positions[7]
		number_missing_5 = list_of_positions[10] - list_of_positions[9]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		list_of_Ns.append("N"*number_missing_3)
		list_of_Ns.append("N"*number_missing_4)
		list_of_Ns.append("N"*number_missing_5)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 14:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		number_missing_3 = list_of_positions[6] - list_of_positions[5]
		number_missing_4 = list_of_positions[8] - list_of_positions[7]
		number_missing_5 = list_of_positions[10] - list_of_positions[9]
		number_missing_6 = list_of_positions[12] - list_of_positions[11]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		list_of_Ns.append("N"*number_missing_3)
		list_of_Ns.append("N"*number_missing_4)
		list_of_Ns.append("N"*number_missing_5)
		list_of_Ns.append("N"*number_missing_6)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)
	if length == 16:
		print(species, length)
		number_missing = list_of_positions[2] - list_of_positions[1]
		number_missing_2 = list_of_positions[4] - list_of_positions[3]
		number_missing_3 = list_of_positions[6] - list_of_positions[5]
		number_missing_4 = list_of_positions[8] - list_of_positions[7]
		number_missing_5 = list_of_positions[10] - list_of_positions[9]
		number_missing_6 = list_of_positions[12] - list_of_positions[11]
		number_missing_7 = list_of_positions[14] - list_of_positions[13]
		list_of_Ns.append("N"*number_missing)
		list_of_Ns.append("N"*number_missing_2)
		list_of_Ns.append("N"*number_missing_3)
		list_of_Ns.append("N"*number_missing_4)
		list_of_Ns.append("N"*number_missing_5)
		list_of_Ns.append("N"*number_missing_6)
		list_of_Ns.append("N"*number_missing_7)
		SpeciesN[species] = list_of_Ns
		print(list_of_Ns)

#Now you need to add the list_of_Ns strings to the end of the sequences. This means we are merging the species_sequence and speciesN dictionaries.
#Use an itertool like zip?

for species in SpeciesN:
	newsequence = [zip(species_sequence[species], SpeciesN[species])]
	newsequenceB = ''.join(str(newsequence))
	species_fused_sequence[species] = newsequenceB
#	print(newsequenceB)
	#why weird annotation? Also, order of strings are wrong.

#Or use a default dictionary? 

#Or write a double loop going through each list of each species in both dictionaries, and adding them to the new dictionary...

#PART 4: PRINTING TO FUSED SEQUENCES TO A NEW FILE

for species in outfile:
	outfile[species] = outfile[species] + "\n" + unicode(str(species_fused_sequence_B[species])) + "\n"
#	print(outfile[species])
	outputFile.write(outfile[species])
