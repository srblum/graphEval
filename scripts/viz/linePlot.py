#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def getNonZeroIndices(inList):
	nonZeroIndices=[]
	for index,value in enumerate(inList):
		if value!=0:
			nonZeroIndices.append(index)
	return nonZeroIndices

def averageNonZero(inList):
	if not all([item==0 for item in inList]):
		total=0
		count=0
		for item in inList:
			if item!=0:
				total+=item
				count+=1
		average=total/count
		return average
	else:
		sys.exit("All zeroes!")


def main():
	colorDict={
		'curoverse':'#a6cee3',
		'debruijn_31k':'#e31a1c',
		'debruijn_63k':'#ff7f00',
		'cactus':'#1f78b4',
		'camel':'#33a02c',
		'simons':'#b2df8a',
		'sbg':'#b15928',
		'prg':'#6a3d9a',
		'refonly':'#000000',
		'trivial':'#b1b300',
		'snp1kg':'#fb9a99',
		'vglr':'#cab2d6',
		'haplo1kg30':'#00ff00',
		'haplo1kg50':'#0000ff',
		'shifted1kg':'#ff0000'
	}
	nameDict={
		'curoverse':'Curoverse',
		'debruijn_31k':'De Bruijn 31',
		'debruijn_63k':'De Bruijn 63',
		'cactus':'Cactus',
		'camel':'Camel',
		'sbg':'7BG',
		'refonly':'Primary',
		'trivial':'Unmerged',
		'vglr':'VGLR',
		'snp1kg':'1KG',
		'haplo1kg30':'1KG Haplo 30',
		'haplo1kg50':'1KG Haplo 50',
		'shifted1kg':'Scrambled',
		'simons':'SGDP',
		'prg':'PRG'
	}
	#Dict of non-N primary reference length for each region
	#Note that there are no N's in the primary reference anyways...
	nonNRefPositionDict={
		'sma':2397625,
		'lrc_kir':1058685,
		'mhc':4970458,
		'brca1':81189,
		'brca2':84989
	}

	def getStats(inFile):
		regionList=['BRCA1','BRCA2','LRC_KIR','MHC','SMA']
		def getZeroList():
			return [0,0,0,0,0]
		statDict=defaultdict(getZeroList)
		algo=''
		with open(inFile) as inFile:
			for line in inFile:
				line=line.strip()
				if line.startswith('$'):
					region=line[1:]
				else:
					if region!='cenx':
						line=line.split('\t')
						#algo,source,positionSum,nonNPositionSum,len(sequenceDict),segmentCount,joinCount,segmentedJoinCount,baseJoinCount
						algo=line[0]
						positions=int(line[2])
						nonNPositions=int(line[3])
						sequences=int(line[4])
						segments=int(line[5])
						joins=float(line[6])
						fullJoins=float(line[7])
						baseJoins=float(line[8])
						if algo in ['shifted1kg','refonly','trivial','prg','curoverse','simons','snp1kg','sbg','cactus','camel','debruijn_63k','vglr','haplo1kg50']:
							# statDict[algo][region]['positions']=positions
							# statDict[algo][region]['nonNPositions']=nonNPositions
							# statDict[algo][region]['joins']=joins
							# statDict[algo][region]['sequences']=sequences
							# statDict[algo][region]['segments']=segments
							# statDict[algo][region]['allJoins']=fullJoins
							# statDict[algo][region]['degree']=round(fullJoins/segments,4)
							# statDict[algo][regionList.index(region.upper())]=round(baseJoins/positions,4)
							# For base degree use baseJoins/positions
							# For compression use nonNRefPositionDict[region]/nonNPositions
							# statDict[algo][regionList.index(region.upper())]=round(nonNRefPositionDict[region]/nonNPositions,4)
							statDict[algo][regionList.index(region.upper())]=nonNRefPositionDict[region]/nonNPositions
		return statDict

	def getCutWidth(inFile):
		dataDict={}
		with open(inFile) as inFile:
			for line in inFile:
				if not line.startswith('\t'):
					line=line.strip().split('\t')
					algo=line[0]
					valueList=[float(x) for x in line[1:]]
					dataDict[algo]=valueList
		return dataDict


	# inFile='cutWidth.tsv'
	# dataDict=getCutWidth(inFile)

	inFile='base_degree_stats_results.txt'
	dataDict=getStats(inFile)
	
	fig, ax = plt.subplots(figsize=(16,6))
	interval=12
	for index,pair in enumerate(sorted(dataDict.items(),reverse=True,key=lambda n:averageNonZero(n[1]))):
		algo,valueList=pair
		nonZeroIndices=np.array(getNonZeroIndices(valueList))
		nonZeroValues=[val for val in valueList if val!=0]
		nonZeroAverage=averageNonZero(valueList)
		plt.plot(nonZeroIndices+4+index*interval,nonZeroValues,color=colorDict[algo],marker='o')	#plot the 5 data points
		plt.plot([4+index*interval,8+index*interval],[nonZeroAverage,nonZeroAverage],color=colorDict[algo]) #plot the average line
		# if algo=='haplo1kg30':
		# 	ax.annotate(nameDict[algo],(5+index*interval-0.4*len(nameDict[algo]),min(nonZeroValues)-0.005),color=colorDict[algo],fontsize=14)
		if index%2==0: offset=-0.1
		else: offset=0
		if min(nonZeroValues)-0.1+offset>0.2:
			ax.annotate(nameDict[algo],(5+index*interval-0.55*len(nameDict[algo]),min(nonZeroValues)-0.1+offset),color=colorDict[algo],fontsize=14)
		else:
			ax.annotate(nameDict[algo],(5+index*interval-0.55*len(nameDict[algo]),max(nonZeroValues)+0.12+offset),color=colorDict[algo],fontsize=14)

	ax.set_ylabel('Compression',fontsize=18)

	#Get rid of the right and top borders
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	plt.tick_params(axis='x',bottom='off',labelbottom='off')

	# plt.show()
	plt.savefig('testCompression.pdf',format='pdf',dpi=400)



if __name__=="__main__":
	main()