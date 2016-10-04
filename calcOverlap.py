## Code to find the overlap between disease gene set and pathway gene sets.
## Written by Sumaiya Nazeen, nazeen@mit.edu
## Python version 2.7.6
## Platform: 3.13.0-35-gneric #62-Ubuntu (64-bit)

## The first input file should be gene sets from MSigDB in gmt format.
## The second input file is the disease gene set selected after applyting BY test on Fisher's combined p-values. It is a tab-separated file with two columns: <gene,p-val>
## The output file should preferrably have a .tab extension.

## Example command to generate the combined results file:
##### $python calcOverlap.py c2.all.v4.0.symbols.gmt sig_asd_by.tab eall_asd.tab 

import sys

def buildPathwayMap(filename):
	lines = [line.strip() for line in open(filename)]
	map = {}
	count = 0
	for line in lines:
		x = line.split('\t')
		if len(x[2:])>=10 and len(x[2:]) <=300: # filter out too small and too large gene sets
			#print x[0], len(x[2:])
			map[x[0].upper()] = x[2:]
			count += 1
	print "Num of selected pathways", count
	return map	

def readGenes(filename):
	lines = [line.strip() for line in open(filename)]
	genes = []
	for line in lines:
		x = line.split('\t')[0]
		y = x.split('///')
		for g in y:
			genes.append(g.upper())
	# print genes
	return genes

def calcEnrichment(map,genes,outfile):
	of = open(outfile,'w')

	for key in map.keys():
		detect = 0
		dlist = []
		for gene in map[key]:
			if gene in genes:
				detect += 1
				dlist.append(gene)
		print key, detect
		if detect > 0:
			s = key + '\t' + str(detect) + '\t' + str(len(map[key])) + '\t' + dlist[0]
			for gene in dlist[1:]:
				s += ','+gene
			s += '\n'
			of.write(s)

	of.close()

def main():
	if len(sys.argv) < 4:
		print "Usage: ./calcOverlap.py <pathway.gmt> <ranklist> <out>"
		sys.exit(1)
	
	pfile = sys.argv[1]
	gfile = sys.argv[2]
	ofile = sys.argv[3]

	map = buildPathwayMap(pfile)
	ranklist = readGenes(gfile)
	calcEnrichment(map,ranklist,ofile)

if __name__ == "__main__":
	main()
