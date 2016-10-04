## Code to combine the output files per pathways collection from enrichment analysis
## Written by Sumaiya Nazeen, nazeen@mit.edu
## Python version 2.7.6
## Platform: 3.13.0-35-gneric #62-Ubuntu (64-bit)

## The output file should preferrably have a .tab extension, so that it can be easily opened with MS Excel for further calculations
## The input file list should contain all the output files for a certain pathway collection specified by the user.

## Example command to generate the combined results file:
##### $python combineOutputOfEnrichment.py kegg/output_combo_kegg_p.tab kegg/phyp_asd.tab kegg/phyp_asthma.tab kegg/phyp_inf.tab kegg/phyp_ckd.tab kegg/phyp_cp.tab kegg/phyp_dc.tab kegg/phyp_ei.tab kegg/phyp_ep.tab kegg/phyp_ibd.tab kegg/phyp_md.tab kegg/phyp_s.tab kegg/phyp_uri.tab

import sys
import numpy as np

def getPathways(inFile):
	lines = [line.strip() for line in open(inFile)]
	paths = []
	for line in lines[1:]:
		paths.append(line.split('\t')[0])
	return paths

def getPs(inFile):
	lines = [line.strip() for line in open(inFile)]
	mF = {}
	for line in lines[1:]:
		x = line.split('\t')
		if x[4] == "NA":
			x[4] = -1
		mF[x[0]] = float(x[4])
	return mF
 
def makeMap(fileList, ofile):
	m = {}
	i = 0
	for f in fileList:
		gF = getPathways(f)
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
		print "Usage: ./combineOutputOfEnrichment.py <outFile> <inFile list>"
		sys.exit(1)

	ofile = sys.argv[1]
	inFiles = sys.argv[2:]

	makeMap(inFiles, ofile)

if __name__ == "__main__":
	main()
