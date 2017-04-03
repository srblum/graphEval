#!/usr/bin/env python

from __future__ import division,print_function
import sys


def getGraphKmerDict(inFile):
	graphKmerDict=dict()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				line=line.strip().split()
				kmer=line[0].upper()
				count=int(line[1])
				graphKmerDict[kmer]=count
	return graphKmerDict

def getReadKmerDict(inFile):
	readKmerDict=dict()
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line and len(line.split())==2:
				line=line.split()
				kmer=line[0]
				count=int(line[1])
				readKmerDict[kmer]=count
	return readKmerDict


def main():
	graphKmerFile,readKmerFile=sys.argv[1:]
	graphKmerDict=getGraphKmerDict(graphKmerFile)
	readKmerDict=getReadKmerDict(readKmerFile)
	# if not all([len(kmer)==20 for kmer in graphKmerDict]):
	# 	print("Not all kmers in graph are length 20!")
	# if not all([len(kmer)==20 for kmer in readKmerDict]):
	# 	print("Not all kmers in reads are length 20!")
	totalReadKmers=sum(readKmerDict.itervalues())
	readKmersInGraph=sum([readKmerDict[read] for read in readKmerDict if read in graphKmerDict])
	print(readKmersInGraph,totalReadKmers)

if __name__=="__main__":
	main()