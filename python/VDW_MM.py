# Created by Dr. Jie Zhang, Manasi Doshi. Copyright 2020 Jie Zhang.
# Updated by Khai Yi Chin

#################################################################################
#				Generate the MM. section in siesta input file					#
# * Update the Grimme Potential C6 and van der Waals radii for the species you  #
# 	need in the following program (SIESTA-4.0 mannual P90). 					#
# * Grimme Potential can be obtained with program "DFT_D3", van der Walls radii	#
#	can be obtained from literatures. (e.g. Semiempirical GGA type density		#
#  	functional constructed with a long-range dispersion correction by Grimme)	#
# * C6AB interaction uses an approximation of: C6AB=(C6AA * C6BB)^0.5 ;  		#
#  	vdW radii for interacting species is the sum of atom vdw radii: Rab=Ra+Rb	#
#################################################################################

# The DFT_D3 program can be obtained from:
# https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/get-the-current-version-of-dft-d3

import sys 
import os
import math
import shutil

######## Update your Date Base here ########
# / atom number / C6AA (eV) / vdW Radii (Bohr)/ #

db=[[6, 1333.083, 0.0321, 'C'],		# C
	[19, 131107.427, 0.0831, 'K'],		# K
	[35, 4548.917, 0.0472, 'Br'],		# Br 
	[79, 9301.375, 0.0605, 'Au']]		# Au
############################################
index_atom=[row[0] for row in db]

print '\tGreeting! Amigo!\nThis code is to generate the vdW force field in Siesta input file.\n'


print '\nPlease input the atom number you are using in SIESTA, follow the index sequence in the block of ChemicalSpeciesLabel (exit with 0):'
index1=[]
counter1=0
while True:             	# Loop continuously
	while True:
		try:
			inp = input("Enter the atom number: ")   	# Get the input
			break
		except Exception:
			print 'Sorry this is not a number, please enter a number'
			pass
	if inp == 0:       	# If it is a blank line...
		break           	# ...break the loop
	if inp not in index_atom:
		print 'Sorry this species is not yet in the data base, please add the potential and vdW radii manually'
		sys.exit(1)
	else:
		index1.append(inp)
		counter1=counter1+1
		print 'The', counter1, 'species is:', db[index_atom.index(inp)][3]

	
print '\nYour index of the atoms will be:', index1,'\n'


ct = open('MM_section','w')
ct.write('MM.UnitsDistance\tBohr\n')
ct.write('MM.UnitsEnergy\tau\n')
ct.write('MM.Grimme.S6\t1.0\n\n')
ct.write('%block MM.Potentials\n')

for x in range(1, counter1+1):
	for y in range(x, counter1+1):
		aa=index_atom.index(index1[x-1]);
		bb=index_atom.index(index1[y-1]);
		cab=math.sqrt(db[aa][1]*db[bb][1])
		vdwr=(db[aa][2]+db[bb][2])
		ct.write("%d %d Grimme %.3f %.3f\n" % (x, y, cab, vdwr))

ct.write('%endblock MM.Potentials\n')

ct.close()



