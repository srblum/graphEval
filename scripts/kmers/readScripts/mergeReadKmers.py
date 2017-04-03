#!/usr/bin/env python
from __future__ import division, print_function
from subprocess import Popen,PIPE
from collections import Counter
from uuid import uuid4
import os, time

def getKmerDict(inFile):
	readKmerDict=Counter()
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line:
				line=line.split()
				kmer=line[-2][-20:]
				count=int(line[-1])
				if count>1:
					readKmerDict[kmer]=count
	return readKmerDict

def mergeDictList(dictList,tempPath):
	processList=[] #List of Popen objects
	inputList=[] #Triplets of [inFile1, inFile2, outFile]
	newDictList=[] #The list of files to be processed next time
	nextInputNo=0

	#Fill inputList with [inFile1, inFile1, outFile] triplets to be processed
	print("Merging list of {} dicts...".format(len(dictList)))
	for i in range(0,len(dictList),4):
		inFiles=dictList[i:i+4]
		if len(inFiles)==4:
			outFile=tempPath+str(uuid4())+'.temp'
			inputList.append([inFiles[0],inFiles[1],inFiles[2],inFiles[3],outFile])
			newDictList.append(outFile)
		else:
			newDictList+=inFiles

	totalInputs=len(inputList)

	processList,nextInputNo=checkRunning(processList,totalInputs,nextInputNo,inputList)
	while len(processList)>0:
		time.sleep(1)
		processList,nextInputNo=checkRunning(processList,totalInputs,nextInputNo,inputList)

	return newDictList


def checkRunning(processList,totalInputs,nextInputNo,inputList,maxProcesses=32):
	"""Check any running processes and start new ones if there are spare slots."""
	for p in range(len(processList)-1,-1,-1): # Check the processes in reverse order
		if processList[p].poll() is not None: # If the process hasn't finished will return None
			del processList[p] # Remove from list - this is why we needed reverse order

	while (len(processList) < maxProcesses) and (nextInputNo < totalInputs): # More to do and some spare slots
		inFile1,inFile2,inFile3,inFile4,outFile=inputList[nextInputNo]
		proc = Popen(['./mergeDict.py',inFile1,inFile2,inFile3,inFile4,outFile])
		processList.append(proc)
		nextInputNo+=1
	return processList, nextInputNo


def writeDict(finalDict,filePath):
	with open(filePath,'w') as outFile:
		for kmer in finalDict:
			outFile.write(kmer+'\t'+str(finalDict[kmer])+'\n')

def main():
	readPath='/hive/users/sblum/ga4gh/reads_k20/'
	tempPath='/scratch/sblum/compareKmers/'

	for region in os.listdir(readPath):
		if region in ['CENX']:
			regionPath=readPath+region+'/'
			outFile=regionPath+'completeThreshold1TotalCounts.dict'
			print("Merging "+region+"...")

			#dictFileList is used to merge together many dicts into one dict
			dictFileList=[regionPath+sample for sample in sorted(os.listdir(regionPath)) if sample.endswith('.jf')]

			#Pairwise merging of dicts until only 2 remain
			while len(dictFileList)>3:
				dictFileList=mergeDictList(dictFileList,tempPath)

			#Merge the last dicts together
			print("Merging the last {} dicts...".format(len(dictFileList)))
			finalDict=Counter()
			for dictFile in dictFileList:
				kmerDict=getKmerDict(dictFile)
				finalDict+=kmerDict
				if dictFile.endswith('.temp'):
					os.remove(dictFile)

			writeDict(finalDict,outFile)


	print("Done.")



if __name__=="__main__":
	main()