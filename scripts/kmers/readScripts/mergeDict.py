#!/usr/bin/env python
from __future__ import division,print_function
import sys, os
from collections import Counter

def getKmerDict(inFile):
	readKmerDict=Counter()
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line:
				line=line.split()
				kmer=line[-2][-20:]
				count=int(line[-1])
				if count>1:
					readKmerDict[kmer]=count
	return readKmerDict


def main():
	inFile1,inFile2,inFile3,inFile4,outFile=sys.argv[1:]
	kmerDict1=getKmerDict(inFile1)
	kmerDict2=getKmerDict(inFile2)
	kmerDict3=getKmerDict(inFile3)
	kmerDict4=getKmerDict(inFile4)
	mergedDict=kmerDict1+kmerDict2+kmerDict3+kmerDict4

	with open(outFile,'w') as outFile:
		for kmer in mergedDict:
			outFile.write(kmer+'\t'+str(mergedDict[kmer])+'\n')

	for inFile in [inFile1, inFile2,inFile3, inFile4]:
		if inFile.endswith('.temp'):
			os.remove(inFile)

if __name__=="__main__":
	main()