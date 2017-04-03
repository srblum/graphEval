#!/usr/bin/env python
from __future__ import division,print_function
import os

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
	revKmerDictFile='/cluster/home/sblum/ga4gh/kmers/graphs/revKmers/LRC_KIR/cactus-lrc_kir_edgeLimit7.dict'
	readFile='/hive/users/sblum/ga4gh/reads/LRC_KIR/HG01494.jf'
	graphKmerDict=getGraphKmerDict(revKmerDictFile)
	totalGraphKmers=sum(graphKmerDict.values())
	totalGraphNKmers=sum([pair[1] for pair in graphKmerDict.items() if 'n' in pair[0]])
	print("Total Kmers in graph:",totalGraphKmers)
	print("Total Kmers in graph with N's in them:",totalGraphNKmers)
	print("Proportion of Kmers in graph with N's in them:",totalGraphNKmers/totalGraphKmers)

	readKmerDict=getReadKmerDict(readFile)
	totalReadKmers=sum(readKmerDict.values())
	totalReadNKmers=sum([pair[1] for pair in readKmerDict.items() if 'N' in pair[0] or 'n' in pair[0]])


	print("Total Kmers in reads:",totalReadKmers)
	print("Total Kmers in reads with N's in them:",totalReadNKmers)
	print("Proportion of Kmers in reads with N's in them:",totalReadNKmers/totalReadKmers)


if __name__=="__main__":
	main()