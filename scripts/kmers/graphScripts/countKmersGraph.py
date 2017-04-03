#!/usr/bin/env python

from __future__ import division, print_function
from subprocess import Popen
from collections import defaultdict
import os, time


def main():
	vgPath='/hive/users/sblum/ga4gh/kmers/graphs/vg/'
	countPath='/hive/users/sblum/ga4gh/kmers/graphs/kmers/'
	maxProcesses=32

	nextInputNo = 0
	inputList=[]
	processList=[]
	global maxProcesses
	global nextInputNo
	global inputList
	global processList

	#Fill inputList with (inFile, outFile) pairs that need to be processed
	for region in ['MHC']:#os.listdir(vgPath):
		inRegionPath=vgPath+region+'/'
		outRegionPath=countPath+region+'/'
		for vgFile in os.listdir(inRegionPath):
			inStem=vgFile.split('.')[0]
			inFile=inRegionPath+vgFile
			outFile=outRegionPath+inStem+'_edgeLimit7.kmer'
			if inStem+'_edgeLimit7.kmer' not in os.listdir(outRegionPath) and inStem.startswith('haplo'):
				inputList.append([inFile,outFile])

	totalInputs = len(inputList)
	global totalInputs

	#NEW CODE
	checkRunning()
	while len(processList)>0:
		time.sleep(1)
		checkRunning()
	print("Done!")



def startNew():
	"""Start a new subprocess if there is work to do """
	global nextInputNo
	global processList

	if nextInputNo < totalInputs:
		inFile,outFile=inputList[nextInputNo]
		proc = Popen(['./countKmersGraph.sh',inFile,outFile])
		print ("Started to Process {}".format(inputList[nextInputNo]))
		nextInputNo+=1
		processList.append(proc)

def checkRunning():
	"""Check any running processes and start new ones if there are spare slots."""
	global processList
	global nextInputNo

	for p in range(len(processList)-1,-1,-1): # Check the processes in reverse order
		if processList[p].poll() is not None: # If the process hasn't finished will return None
			del processList[p] # Remove from list - this is why we needed reverse order

	while (len(processList) < maxProcesses) and (nextInputNo < totalInputs): # More to do and some spare slots
		startNew()



if __name__=="__main__":
	main()