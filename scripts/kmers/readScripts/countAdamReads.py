#!/usr/bin/env python
from __future__ import division, print_function
from subprocess import Popen
import os,time




def maximizeProcesses(processList,totalInputs,nextInputNo,inputList,maxProcesses):
	"""Check any running processes and start new ones if there are spare slots."""
	for p in range(len(processList)-1,-1,-1): # Check the processes in reverse order
		if processList[p].poll() is not None: # If the process hasn't finished will return None
			del processList[p] # Remove from list - this is why we needed reverse order

	while (len(processList) < maxProcesses) and (nextInputNo < totalInputs): # More to do and some spare slots
		inFile,outFile=inputList[nextInputNo]
		proc = Popen(['./countAdamReads.sh',inFile,outFile])
		processList.append(proc)
		nextInputNo+=1
	return processList,nextInputNo


def main():
	#maxProcesses dictates the number of simultaneous processes that will be run
	maxProcesses=32
	inReadPath='/hive/users/sblum/ga4gh/reads/'
	outReadPath='/hive/users/sblum/ga4gh/reads_k20/'
	for region in os.listdir(inReadPath):
		if region=='SMA':
			print("Processing "+region+"...")
			inRegionPath=inReadPath+region+'/'
			outRegionPath=outReadPath+region+'/'

			#Make list of [inFile, outFile] pairs that have not already been processed
			inputList=[]

			for inFile in os.listdir(inRegionPath):
				if inFile.endswith('.rev.fa'):
					sample=inFile.split('.')[0]
					if sample+'.jf' not in os.listdir(outRegionPath):
						inFile=inRegionPath+inFile
						outFile=outRegionPath+sample+'.jf'
						inputList.append([inFile,outFile])

			totalInputs=len(inputList)
			nextInputNo=0

			processList=[]

			processList,nextInputNo=maximizeProcesses(processList,totalInputs,nextInputNo,inputList,maxProcesses)
			while len(processList)>0:
				time.sleep(1)
				processList,nextInputNo=maximizeProcesses(processList,totalInputs,nextInputNo,inputList,maxProcesses)
			print("Done!")



if __name__=="__main__":
	main()