#!/usr/bin/env python

from __future__ import division, print_function
from subprocess import Popen
import os

def main():
	#This number dictates the number of simultaneous processes that will be run
	#Should not exceed the number of cores on your machine
	numCores=32
	processList=[]
	# This is new path.  Old path was '/hive/users/anovak/ga4gh/bake-off/hgvm/reads_all/'
	adamReadPath='/cluster/home/anovak/hive/ga4gh/bake-off/hgvm/reads_all/'
	seanReadPath='/hive/users/sblum/ga4gh/reads/'
	for region in os.listdir(adamReadPath):
		if region != 'good.txt':
			print("Processing "+region+"...")
			regionPath=adamReadPath+region+'/'
			fastaPath=seanReadPath+region+'/'
			count=0
			for index in range(0,len(os.listdir(regionPath)),numCores):
				for sampleDir in sorted(os.listdir(regionPath))[index:index+numCores]:
					fastqPath=regionPath+sampleDir+'/'
					if sampleDir+'.fa' not in os.listdir(fastaPath):
						count+=1
						command=['./getAdamReads.sh',fastqPath+sampleDir+'.bam.fq',fastaPath+sampleDir+'.fa']
						processList.append(Popen(command))
				exitCodes=[p.wait() for p in processList]
			print("Transferred {} files...".format(count))
	print('Finished.')

if __name__=="__main__":
	main()