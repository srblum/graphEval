#!/usr/bin/env python

from __future__ import division,print_function
import sys

def getKmerSet(inFile):
	kmerSet=set()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				kmerSet.add(line.strip())
	return kmerSet


def main():
	inFile1,inFile2,outFile=sys.argv[1:]
	kmerSet1=getKmerSet(inFile1)
	kmerSet2=getKmerSet(inFile2)
	with open(outFile,'w') as outFile:
		for kmer in kmerSet1|kmerSet2:
			outFile.write(kmer+'\n')

if __name__=="__main__":
	main()