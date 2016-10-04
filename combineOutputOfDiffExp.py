## Code to combine the output files per disease from individual differential expression analysis
## Written by Sumaiya Nazeen, nazeen@mit.edu
## Python version 2.7.6
## Platform: 3.13.0-35-gneric #62-Ubuntu (64-bit)

## The output file should preferrably have a .tab extension, so that it can be easily opened with MS Excel for further calculations
## The input file list should contain all the output files under a certain disease specified by the user.

## Example command to generate the combined results file:
##### $python combineOutputOfDiffExp.py combo_asthma.tab diff_a1.tab diff_a2.tab diff_a3.tab diff_a4.tab diff_a5.tab

import sys
import numpy as np

def getGenes(inFile):
	lines = [line.strip() for line in open(inFile)]
	genes = []
	for line in lines[1:]:
		genes.append(line.split('\t')[1])
	return genes

def getPs(inFile):
	lines = [line.strip() for line in open(inFile)]
	mF = {}
	for line in lines[1:]:
		x = line.split('\t')
		if x[3] == "NA":
			x[3] = -1
		mF[x[1]] = float(x[3])
	return mF
 
def makeMap(fileList, ofile):
	m = {}
	i = 0
	for f in fileList:
		gF = getGenes(f)
		for g in gF:
			if g not in m.keys(): 	
				m[g] = [-1]*len(fileList)
	
		mF = getPs(f)
		for g in gF:
			m[g][i] = mF[g]
	
		i = i+1
	
	of = open(ofile,'w')

	for g in m.keys():
		s = g
		for j in xrange(len(m[g])):
			s += '\t' + str(m[g][j])
		s += '\n'
		of.write(s)
	
	of.close()

def main():
	if len(sys.argv) < 3:
		print "Usage: ./combineOutputOfDiffExp.py <outFile> <inFile list>"
		sys.exit(1)

	ofile = sys.argv[1]
	inFiles = sys.argv[2:]

	makeMap(inFiles, ofile)

if __name__ == "__main__":
	main()
