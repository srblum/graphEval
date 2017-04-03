#!/usr/bin/env python

from __future__ import division,print_function
import sys


def getKmerDict(inFile):
	kmerDict=dict()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				line=line.strip().split()
				kmer=line[0]
				count=int(line[1])
				kmerDict[kmer]=count
	return kmerDict



def main():
	graphKmerFile,readKmerFile,region=sys.argv[1:]
	graphKmerDict=getKmerDict(graphKmerFile)
	readKmerDict=getKmerDict(readKmerFile)
	totalReadKmerInstances=sum(readKmerDict.itervalues())
	totalGraphKmerInstances=sum(graphKmerDict.itervalues())

	readKmerInstancesInGraph=0
	graphKmerInstancesInReads=0
	uniqueReadKmersInstancesInGraph=0
	for kmer in readKmerDict:
		if kmer in graphKmerDict:
			readKmerInstancesInGraph+=readKmerDict[kmer]
			graphKmerInstancesInReads+=graphKmerDict[kmer]
			if graphKmerDict[kmer]==1:
				uniqueReadKmersInstancesInGraph+=readKmerDict[kmer]

	instancePrecision=graphKmerInstancesInReads/totalGraphKmerInstances
	instanceRecall=readKmerInstancesInGraph/totalReadKmerInstances
	uniqueRecall=uniqueReadKmersInstancesInGraph/totalReadKmerInstances
	print(instancePrecision,instanceRecall,uniqueRecall,len(graphKmerDict),graphKmerFile,region)

if __name__=="__main__":
	main()