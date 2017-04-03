#!/usr/bin/env python

from __future__ import division,print_function
import sys


def getGraphKmerSet(inFile):
	graphKmerSet=set()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				line=line.strip().split()
				kmer=line[0].upper()
				graphKmerSet.add(kmer)
	return graphKmerSet

def getReadKmerSet(inFile):
	readKmerSet=set()
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line and len(line.split())==2:
				line=line.split()
				kmer=line[0]
				readKmerSet.add(kmer)
	return readKmerSet


def main():
	args=sys.argv
	graphKmerFile,readKmerFile,outFile=args[1:4]
	graphKmerSet=getGraphKmerSet(graphKmerFile)
	readKmerSet=getReadKmerSet(readKmerFile)
	with open(outFile,'w') as outFile:
		for kmer in graphKmerSet&readKmerSet:
			outFile.write(kmer+'\n')

	#If a fourth argument was passed, write all the unfound kmers
	#To a file named using the fourth argument
	if len(args)==5:
		newGraphKmerFile=args[4]
		with open(newGraphKmerFile,'w') as outFile:
			for kmer in graphKmerSet:
				if kmer not in readKmerSet:
					outFile.write(kmer+'\n')

if __name__=="__main__":
	main()