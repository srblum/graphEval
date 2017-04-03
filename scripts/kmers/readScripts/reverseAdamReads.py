#!/usr/bin/env python

from __future__ import division, print_function
from subprocess import Popen
import os

def main():
	#This number dictates the number of simultaneous processes that will be run
	#Should not exceed the number of cores on your machine
	numCores=32
	processList=[]
	seanReadPath='/hive/users/sblum/ga4gh/reads/'
	for region in os.listdir(seanReadPath):
		print("Processing "+region+"...")
		regionPath=seanReadPath+region+'/'
		for index in range(0,len(os.listdir(regionPath)),numCores):
			for inFile in sorted(os.listdir(regionPath))[index:index+numCores]:
				if not inFile.endswith('rev.fa') and not inFile.endswith('.jf'):
					sample=inFile.split('.')[0]
					outFile=regionPath+sample+'.rev.fa'
					inFile=regionPath+inFile
					if sample+'rev.fa' not in os.listdir(regionPath):
						command=['./reverseAdamReads.sh',inFile,outFile]
						processList.append(Popen(command))
			exitCodes=[p.wait() for p in processList]
	print('Finished.')

if __name__=="__main__":
	main()