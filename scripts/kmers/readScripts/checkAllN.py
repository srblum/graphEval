#!/usr/bin/env python
from __future__ import division,print_function
from subprocess import Popen
import os


def main():
	#This number dictates the number of simultaneous processes that will be run
	#Should not exceed the number of cores on your machine
	numCores=32
	seanReadPath='/hive/users/sblum/ga4gh/reads/'
	for region in os.listdir(seanReadPath):
		print("Processing "+region+"...")
		regionPath=seanReadPath+region+'/'

		#Make list of .jf files to check for N's in
		sampleList=[]
		for inFile in os.listdir(regionPath):
			if inFile.endswith('.jf'):
				sample=inFile.split('.')[0]
				sampleList.append(sample)


		#Process files
		for index in range(0,len(sampleList),numCores):
			processList=[]
			for sample in sampleList[index:index+numCores]:
					inFile=regionPath+sample+'.jf'
					command=['./checkN.py',inFile]
					processList.append(Popen(command))
			exitCodes=[p.wait() for p in processList]
if __name__=="__main__":
	main()