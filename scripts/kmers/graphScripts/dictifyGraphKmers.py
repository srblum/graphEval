#!/usr/bin/env python
from __future__ import division,print_function
import os
from subprocess import Popen

def main():
	numCores=32
	graphPath='/hive/users/sblum/ga4gh/kmers/graphs/revKmers/'

	inputList=[]
	for region in ['BRCA1','BRCA2','SMA','MHC','LRC_KIR']:
		regionPath=graphPath+region+'/'
		for kmerFile in os.listdir(regionPath):
				inStem=kmerFile.split('.')[0]
				inFile=regionPath+kmerFile
				outFile=regionPath+inStem+'.dict'
				if inStem+'.dict' not in os.listdir(regionPath) and inStem.endswith('edgeLimit7'):
					inputList.append([inFile,outFile])

	for index in range(0,len(inputList),numCores):
		print("Processing reverse kmer files {} through {}...".format(index+1,min([index+numCores,len(inputList)])))
		processList=[]
		for inFile,outFile in inputList[index:index+numCores]:
			command=['./dictify.py',inFile,outFile]
			processList.append(Popen(command))
		exitcodes=[p.wait() for p in processList]
	print("Done.")


if __name__=="__main__":
	main()
