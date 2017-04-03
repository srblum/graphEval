#!/usr/bin/env python

from __future__ import division, print_function
from subprocess import Popen
from collections import defaultdict
import os

def getServers(inFile):
	serverDict=defaultdict(list)
	with open(inFile) as inFile:
		for line in inFile:
			if not line.startswith('#'):
				line=line.strip().split('\t')
				region,url,algo,source=line[:4]
				if region!='region':
					region=region.upper()
					serverDict[region].append([url,algo,source])
	return serverDict

def main():
	numCores=32
	serverFile='/cluster/home/sblum/ga4gh/hgvm-graph-bakeoff-evalutations/graph_servers_azure.tsv'
	serverDict=getServers(serverFile)
	processList=[]
	for region in serverDict:
		print("Processing "+region+"...")
		graphPath='/hive/users/sblum/ga4gh/kmers/graphs/vg/'+region+'/'
		for url,algo,source in serverDict[region]:
			graph=url.split('/')[-2]
			url+='v0.6.g'
			outFile=graphPath+graph+'.vg'
			command=['./getGraphs.sh',url,outFile]
			if graph+'.vg' not in os.listdir(graphPath):
				command=['./getGraphs.sh',url,outFile]
				processList.append(Popen(command))
				if len(processList)==numCores:
					exitCodes=[p.wait() for p in processList]
					processList=[]
	if processList:
		exitCodes=[p.wait() for p in processList]
	print('Finished.')

if __name__=="__main__":
	main()