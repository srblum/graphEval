#!/usr/bin/env python
from __future__ import division, print_function
from subprocess import Popen,PIPE
from collections import defaultdict
import os, time


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

def maximizeProcesses(processList,totalInputs,nextInputNo,inputList,resultDict,maxProcesses=32):
	"""Check any running processes and start new ones if there are spare slots."""
	for p in range(len(processList)-1,-1,-1): # Check the processes in reverse order
		process=processList[p]
		if process.poll() is not None: # If the process hasn't finished will return None
			result=process.communicate()[0]
			instancePrecision,instanceRecall,uniqueRecall,graphKmerNo,graphKmerFile,region=result.strip().split()
			print("There are {} kmers in {}".format(graphKmerNo,graphKmerFile))
			fileStem='_'.join(graphKmerFile.split('/')[-1].split('_')[:-1])
			resultDict[region][fileStem]=[instancePrecision,instanceRecall,uniqueRecall]
			del processList[p] # Remove from list - this is why we needed reverse order

	while (len(processList) < list) and (nextInputNo < totalInputs): # More to do and some spare slots
		graphKmerFile,readKmerFile,region=inputList[nextInputNo]
		print("Processing "+graphKmerFile+"...")
		proc = Popen(['./getOverlap.py',graphKmerFile,readKmerFile,region],stdout=PIPE,stdin=PIPE)
		processList.append(proc)
		nextInputNo+=1
	return processList, nextInputNo, resultDict

def main():
	readPath='/hive/users/sblum/ga4gh/reads_k20/'
	graphPath='/hive/users/sblum/ga4gh/kmers/graphs/revKmers/'
	resultFile='/hive/users/sblum/ga4gh/kmers/oct08_haplo50_Unique_Threshold1Results.txt'
	readFileName='completeThreshold1TotalCounts.dict'

	resultDict=defaultdict(dict)

	processList=[]
	inputList=[]

	#Fill inputList with (graphKmerFile, readKmerFile) pairs
	for region in os.listdir(graphPath):
		graphRegionPath=graphPath+region+'/'
		readRegionPath=readPath+region+'/'
		readKmerFile=readRegionPath+readFileName
		for graphKmerFile in os.listdir(graphRegionPath):
			if ('haplo1kg50' in graphKmerFile) and graphKmerFile.endswith('edgeLimit7.dict'):
				graphKmerFile=graphRegionPath+graphKmerFile
				inputList.append((graphKmerFile,readKmerFile,region))

	totalInputs=len(inputList)
	nextInputNo=0

	processList,nextInputNo,resultDict=maximizeProcesses(processList,totalInputs,nextInputNo,inputList,resultDict)
	while len(processList)>0:
		processList,nextInputNo,resultDict=maximizeProcesses(processList,totalInputs,nextInputNo,inputList,resultDict)
		time.sleep(1)


	with open(resultFile,'w') as outFile:
		for region in resultDict:
			outFile.write('$'+region+'\n')
			for graph in resultDict[region]:
				instancePrecision,instanceRecall,uniqueRecall=resultDict[region][graph]
				outFile.write(graph+'\t'+instancePrecision+'\t'+instanceRecall+'\t'+uniqueRecall+'\n')
	print("Done.")

if __name__=="__main__":
	main()