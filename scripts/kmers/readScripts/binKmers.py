#!/usr/bin/env python
from __future__ import print_function, division
import os
from collections import defaultdict

def getCounts(inFile):
	countDict=defaultdict(int)
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip().split()
			count=int(line[1])
			countDict[count]+=1
	return countDict

def main():
	readPath='/hive/users/sblum/ga4gh/reads_k20/'
	outFile='/cluster/home/sblum/ga4gh/kmers/threshold1ReadKmerBins.txt'
	with open(outFile,'w') as outFile:
		for region in os.listdir(readPath):
			print ("Processing "+region+"...")
			outFile.write('@'+region+'\n')
			regionPath=readPath+region+'/'
			dictFile=regionPath+'threshold1TotalCounts.dict'
			countDict=getCounts(dictFile)
			for count in sorted(countDict):
				outFile.write(str(count)+'\t'+str(countDict[count])+'\n')


if __name__ == "__main__":
	main()