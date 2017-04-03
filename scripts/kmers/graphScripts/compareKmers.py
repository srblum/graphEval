#!/usr/bin/env python
from __future__ import division, print_function
from subprocess import Popen,PIPE
from collections import defaultdict
from uuid import uuid4
import os


def getKmerSet(inFile):
	kmerSet=set()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				kmerSet.add(line.strip())
	return kmerSet

def getGraphKmerDict(inFile):
	graphKmerDict=dict()
	with open(inFile) as inFile:
		for line in inFile:
			if line:
				line=line.strip().split()
				kmer=line[0].upper()
				count=int(line[1])
				graphKmerDict[kmer]=count
	return graphKmerDict

def mergeSetList(setList,tempPath):
	processList=[]
	newSetList=[]
	print("Merging list of {} sets...".format(len(setList)))
	for i in range(0,len(setList),2):
		inFiles=setList[i:i+2]
		if len(inFiles)==2:
			outFile=tempPath+str(uuid4())+'.set'
			command=['./mergeSet.py',inFiles[0],inFiles[1],outFile]
			processList.append(Popen(command,stdout=PIPE,stdin=PIPE))
			newSetList.append(outFile)
		else:
			newSetList.append(inFiles[0])
		if len(processList)==32:
			[p.wait() for p in processList]
			processList=[]
	[p.wait() for p in processList]
	return newSetList


def main():
	readPath='/hive/users/sblum/ga4gh/reads/'
	graphPath='/cluster/home/sblum/ga4gh/kmers/graphs/revKmers/'
	tempPath='/scratch/sblum/compareKmers/'
	resultFile='/cluster/home/sblum/ga4gh/kmers/LRC_KIR_results.txt'

	def nestedDD():
		return defaultdict(dict)

	resultDict=defaultdict(nestedDD)

	readSubSample=False

	#Compute recall
	for region in os.listdir(graphPath):
		if region=='LRC_KIR':
			print("Processing recall for "+region)
			graphRegionPath=graphPath+region+'/'
			readRegionPath=readPath+region+'/'
			for graphKmerFile in os.listdir(graphRegionPath):
				if graphKmerFile.endswith('.dict'):
					print(graphKmerFile)
					overlapSum=0
					totalSum=0
					sampleList=[sample for sample in sorted(os.listdir(readRegionPath)) if sample.endswith('.jf')]
					if readSubSample:
						sampleList=sampleList[:readSubSample]
					for i in range(0,len(sampleList),32):
						print("Processing recall for samples {} through {}...".format(i+1,i+32))
						processList=[]
						for readKmerFile in sampleList[i:i+32]:
							graphFile=graphRegionPath+graphKmerFile
							readFile=readRegionPath+readKmerFile
							command=['./compareRecall.py',graphFile,readFile]
							processList.append(Popen(command,stdout=PIPE,stdin=PIPE))
						resultList=[p.communicate()[0] for p in processList]
						
						for result in resultList:
							result=result.strip().split()
							overlap,total=result
							overlapSum+=int(overlap)
							totalSum+=int(total)
					recall=overlapSum/totalSum
					print("Recall: "+str(recall))
					fileStem='_'.join(graphKmerFile.split('_')[:-1])
					resultDict[region][fileStem]['recall']=recall

	#Compute Precision
	for region in os.listdir(graphPath):
		if region=='LRC_KIR':
			print("Processing precision for "+region)
			graphRegionPath=graphPath+region+'/'
			readRegionPath=readPath+region+'/'
			for graphKmerFile in os.listdir(graphRegionPath):
				if graphKmerFile.endswith('.dict'):
					print(graphKmerFile)
					fileStem='_'.join(graphKmerFile.split('_')[:-1])
					sampleList=[sample for sample in sorted(os.listdir(readRegionPath)) if sample.endswith('.jf')]
					if readSubSample:
						sampleList=sampleList[:readSubSample]

					#Map kmers in graph to counts to compute precision in the end
					graphKmerDict=getGraphKmerDict(graphRegionPath+graphKmerFile)

					#setFileList is just a list of files to delete at the end of each graph
					setFileList=[]

					#tempSetFileList is used to merge together many sets into one set
					tempSetFileList=[]

					#A counter to track which batch of samples we're processing
					batchCount=0
					for i in range(0,len(sampleList),32):
						print("Processing precision for samples {} through {}...".format(i+1,i+32))

						#processList will contain 32 commands at a time, each one calling
						#comparePrecision.py on a graph and a sample
						processList=[]

						#newGraphKmerFile will contain only the kmers NOT found in 
						#previous batches and add this file to the list of files to delete at the end
						newGraphKmerFile=tempPath+str(batchCount)+graphKmerFile
						setFileList.append(newGraphKmerFile)

						#If we're on the first batch, use the original graph kmer file for comparison
						#Otherwise use the new graph kmer file from the previous batch
						if i==0:
							graphFile=graphRegionPath+graphKmerFile
						else:
							graphFile=tempPath+str(batchCount-1)+graphKmerFile
						batchCount+=1

						#Only fill a newGraphKmerFile once, then set this flag to False
						trackUnfoundKmers=True
						for readKmerFile in sampleList[i:i+32]:
							readFile=readRegionPath+readKmerFile

							#setFile contains the set of all kmers shared by the graph and sample
							setFile=tempPath+fileStem+readKmerFile+'.set'
							setFileList.append(setFile)
							tempSetFileList.append(setFile)

							#Use the first sample in each batch of 32 to generate a new graph kmer file
							#to use in subsequent calls to comparePrecision.py
							if trackUnfoundKmers:
								command=['./comparePrecision.py',graphFile,readFile,setFile,newGraphKmerFile]
								trackUnfoundKmers=False
							else:
								command=['./comparePrecision.py',graphFile,readFile,setFile]

							processList.append(Popen(command))
						[p.wait() for p in processList]
					print("Merging sets...")

					#Pairwise merging of sets until only 2 remain
					while len(tempSetFileList)>2:
						tempSetFileList=mergeSetList(tempSetFileList,tempPath)
						setFileList+=tempSetFileList

					#Merge the last two sets together
					finalSet=getKmerSet(tempSetFileList[0])|getKmerSet(tempSetFileList[1])

					for setFile in setFileList:
						if os.path.isfile(setFile):
							os.remove(setFile)

					totalOverlap=sum([graphKmerDict[kmer] for kmer in graphKmerDict if kmer in finalSet])
					totalSum=sum(graphKmerDict.itervalues())
					precision=totalOverlap/totalSum
					print("Precision: "+str(precision))
					resultDict[region][fileStem]['precision']=precision

	with open(resultFile,'w') as outFile:
		for region in resultDict:
			outFile.write('$'+region+'\n')
			for graph in resultDict[region]:
				recall=resultDict[region][graph]['recall']
				precision=resultDict[region][graph]['precision']
				outFile.write(graph+'\t'+str(precision)+'\t'+str(recall)+'\n')

if __name__=="__main__":
	main()