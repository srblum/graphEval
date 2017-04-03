#!/usr/bin/env python
from __future__ import division, print_function
from subprocess import Popen
from collections import defaultdict
import os

def main():
	numCores=32
	kmerPath='/hive/users/sblum/ga4gh/kmers/graphs/kmers/'
	revPath='/hive/users/sblum/ga4gh/kmers/graphs/revKmers/'
	inputList=[]
	for region in ['BRCA1','BRCA2','SMA','LRC_KIR','MHC']:
		inRegionPath=kmerPath+region+'/'
		outRegionPath=revPath+region+'/'
		for kmerFile in os.listdir(inRegionPath):
			inStem=kmerFile.split('.')[0]
			inFile=inRegionPath+kmerFile
			outFile=outRegionPath+inStem+'.rev'
			if inStem+'.rev' not in os.listdir(outRegionPath) and inStem.endswith('edgeLimit7'):
				inputList.append([inFile,outFile])

	for index in range(0,len(inputList),numCores):
		print("Processing kmer files {} through {}...".format(index+1,min([index+32,len(inputList)])))
		processList=[]
		for inFile,outFile in inputList[index:index+numCores]:
			command=['./reverseComp.py','--graph','--input',inFile,'--output',outFile]
			processList.append(Popen(command))
		exitCodes=[p.wait() for p in processList]
	print('Finished.')

if __name__=="__main__":
	main()