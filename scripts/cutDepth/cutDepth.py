#!/usr/bin/env python
from __future__ import division, print_function
import json, os
from collections import defaultdict


def main():
	inPath='/cluster/home/sblum/ga4gh/cutDepth/sortedVG/'
	outPath='/cluster/home/sblum/ga4gh/cutDepth/results/'
	summaryOutFile='/cluster/home/sblum/ga4gh/cutDepth/results.txt'

	for region in os.listdir(inPath):
		if region!='CENX':

			inRegionPath=inPath+region+'/'
			outRegionPath=outPath+region+'/'
			for inFile in os.listdir(inRegionPath):

				#Chop off the trailing "-sorted.json" from the inFile name
				baseFile=inFile[:-12]
				sampleOutFile=baseFile+'.tsv'

				if 'haplo1kg50' in baseFile or 'snp1kg' in baseFile:

					print("Processing "+baseFile+"...")

					with open(outRegionPath+sampleOutFile,'w') as sampleOutFile:
						sampleOutFile.write("@"+baseFile+"\t")
						with open(inRegionPath+inFile) as inFile:
							text=inFile.read()
							#A dict with keys 'edge', 'node', and 'path'
							#Each node has 'id', 'name', and 'sequence'
							#Each edge has 'to', 'from', and maybe 'to_end' and 'from_start'
							json_=json.loads(text)

							#Dictionary mapping node ID's to lengths
							lengthDict={}
							if 'edge' in json_ and 'node' in json_:

								#Get the lengths of all nodes and store them in lengthDict
								#This is important so we can factor in the lengths of nodes
								#When computing cut width
								for node in json_['node']:
									id_=int(node['id'])
									lengthDict[id_]=len(node['sequence'])

								gapSum=sum(lengthDict.values())-1

								#Special edge type breakdown
								#Not important for general cut width
								standardLengthDict=defaultdict(int)
								feedbackLengthDict=defaultdict(int)
								reverseLengthDict=defaultdict(int)
								standardEdgeSum=0
								reverseEdgeSum=0
								feedbackEdgeSum=0
								standardCount=0
								feedbackCount=0
								reverseCount=0


								for edge in json_['edge']:
									source=int(edge['from'])
									target=int(edge['to'])
									backward=True if source>target else False
									toEnd=True if 'to_end' in edge else False
									fromStart=True if 'from_start' in edge else False

									#Old length
									# length=source-target if backward else target-source

									#New length
									if source>target:
										backward=True
										length=1
										for id_ in range(target+1,source):
											length+=lengthDict[id_]
									elif source==target:
										if not toEnd and not fromStart:
											length=lengthDict[source]-1
											backward=True
										else:
											backward=False
											length=0
									else:
										backward=False
										length=1
										for id_ in range(source+1,target):
											length+=lengthDict[id_]


									logLength=length if length<2 else len(str(length))
									#Reversing edge
									if (fromStart and not toEnd) or (toEnd and not fromStart):
										reverseEdgeSum+=length
										reverseCount+=1
										reverseLengthDict[logLength]+=1
									#Feedback edge
									elif fromStart and toEnd and not backward:
										feedbackEdgeSum+=length
										feedbackCount+=1
										feedbackLengthDict[logLength]+=1
									elif not fromStart and not toEnd and backward:
										feedbackEdgeSum+=length
										feedbackCount+=1
										feedbackLengthDict[logLength]+=1
									#Standard edge
									elif fromStart and toEnd and backward:
										standardEdgeSum+=length
										standardCount+=1
										standardLengthDict[logLength]+=1
									elif not fromStart and not toEnd and not backward:
										standardEdgeSum+=length
										standardCount+=1
										standardLengthDict[logLength]+=1
									else:
										print ('Uncaught edge.')

								baseEdges=sum([length-1 for length in lengthDict.itervalues()])
								edgeSum=standardEdgeSum+reverseEdgeSum+feedbackEdgeSum+baseEdges
								standardPortion=round(standardEdgeSum/edgeSum,2)
								feedbackPortion=round(feedbackEdgeSum/edgeSum,2)
								reversePortion=round(reverseEdgeSum/edgeSum,2)
								averageCutWidth=round(edgeSum/gapSum,20)
								summaryString=str(averageCutWidth)+'\t'+str(edgeSum)+'\t'+str(standardPortion)+'\t'+str(feedbackPortion)+'\t'+str(reversePortion)+'\t'+str(standardCount)+'\t'+str(feedbackCount)+'\t'+str(reverseCount)+'\n'
								sampleOutFile.write(summaryString)
								sampleOutFile.write('@Standard\n')
								for length,count in sorted(standardLengthDict.items(),key=lambda n:n[0]):
									sampleOutFile.write(str(length)+'\t'+str(count)+'\n')
								sampleOutFile.write('@Feedback\n')
								for length,count in sorted(feedbackLengthDict.items(),key=lambda n:n[0]):
									sampleOutFile.write(str(length)+'\t'+str(count)+'\n')
								sampleOutFile.write('@Reverse\n')
								for length,count in sorted(reverseLengthDict.items(),key=lambda n:n[0]):
									sampleOutFile.write(str(length)+'\t'+str(count)+'\n')
							else:
								sampleOutFile.write("NULL")

if __name__=="__main__":
	main()