#!/usr/bin/env python

from __future__ import division,print_function
import sys
from collections import defaultdict

def getDict(inFile):
	kmerDict=defaultdict(int)
	with open(inFile) as inFile:
		for line in inFile:
			if not line.startswith('>'):
				if line and not 'N' in line and not 'n' in line:
					line=line.strip().upper()
					kmerDict[line]+=1
	return kmerDict


def writeDict(kmerDict, outFile):
	with open(outFile,'w') as outFile:
		for kmer in kmerDict:
			outFile.write(kmer+'\t'+str(kmerDict[kmer])+'\n')

args=sys.argv
inFile,outFile=args[1:3]
kmerDict=getDict(inFile)
writeDict(kmerDict,outFile)

