#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt


# def normalize(statDict,regionList,statList):
# 	"""
# 	Returns a dict of 
# 	algo: [averageNormalizedSequences,averageNormalizedSegments,averageNormalizedJoins,averageNormalizedFullJoins]
# 	"""
# 	def getDictDict():
# 		return defaultdict(dict)
# 	normalizedStatDict=defaultdict(getDictDict)
# 	maxDict=defaultdict(dict)
# 	averageNormalizedStatDict=defaultdict(dict)
# 	for stat in statList:
# 		for region in regionList:
# 			algoList=[]
# 			for algo in statDict:
# 				if region in statDict[algo]:
# 					algoList.append(statDict[algo][region][stat])
# 			maxStat=max(algoList)
# 			maxDict[region][stat]=maxStat
# 	for algo in statDict:
# 		for stat in statList:
# 			for region in statDict[algo]:
# 				normalizedStatDict[algo][region][stat]=statDict[algo][region][stat]/maxDict[region][stat]
# 	for algo in normalizedStatDict:
# 		for stat in statList:
# 			total=0
# 			for region in normalizedStatDict[algo]:
# 				total+=normalizedStatDict[algo][region][stat]
# 			average=total/len(normalizedStatDict[algo])
# 			averageNormalizedStatDict[algo][stat]=average
# 	return averageNormalizedStatDict

def main():
	colorDict={
		'curoverse':'#a6cee3',
		'debruijn_k31':'#e31a1c',
		'debruijn_k63':'#ff7f00',
		'cactus':'#1f78b4',
		'camel':'#33a02c',
		'simons':'#b2df8a',
		'sbg':'#b15928',
		'prg':'#6a3d9a',
		'refonly':'#000000',
		'trivial':'#ffff99',
		'snp1000g':'#fb9a99',
		'vglr':'#cab2d6',
		'vg_haplo1kg':'#fdbf6f'
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
		'shifted1kg':'Control',
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
		def getIntDict():
			return defaultdict(int)
		def getDictDict():
			return defaultdict(getIntDict)
		statDict=defaultdict(getDictDict)
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
						if algo in ['prg','snp1kg','sbg','cactus','camel','debruijn_63k','haplo1kg50','shifted1kg']:
							statDict[algo][region]['positions']=positions
							statDict[algo][region]['nonNPositions']=nonNPositions
							statDict[algo][region]['joins']=joins
							statDict[algo][region]['sequences']=sequences
							statDict[algo][region]['segments']=segments
							statDict[algo][region]['allJoins']=fullJoins
							statDict[algo][region]['degree']=round(fullJoins/segments,4)
							statDict[algo][region]['baseDegree']=round(baseJoins/positions,4)
							statDict[algo][region]['compression']=round(nonNRefPositionDict[region]/nonNPositions,4)
		return statDict


		
	statList=['positions','nonNPositions','sequences','segments','joins','allJoins','degree','baseDegree','compression']
	regionList=['brca1','brca2','lrc_kir','mhc','sma']
	inFile="base_degree_stats_results.txt"
	statDict=getStats(inFile)
	# averageNormalizedStatDict=normalize(statDict,regionList,statList)

	# N=5
	# ind = np.arange(N)*2+0.2  # the x locations for the groups
	# width = 0.15       # the width of the bars

	# fig, ax = plt.subplots()

	# cactusRects=ax.bar(ind+width,[averageNormalizedStatDict['cactus'][stat] for stat in statList],width,color=colorDict['cactus'])
	# camelRects=ax.bar(ind+2*width,[averageNormalizedStatDict['camel'][stat] for stat in statList],width,color=colorDict['camel'])
	# k31Rects=ax.bar(ind+3*width,[averageNormalizedStatDict['debruijn_31k'][stat] for stat in statList],width,color=colorDict['debruijn_k31'])
	# k63Rects=ax.bar(ind+4*width,[averageNormalizedStatDict['debruijn_63k'][stat] for stat in statList],width,color=colorDict['debruijn_k63'])
	# prgRects=ax.bar(ind+5*width,[averageNormalizedStatDict['prg'][stat] for stat in statList],width,color=colorDict['prg'])
	# sbgRects=ax.bar(ind+6*width,[averageNormalizedStatDict['sbg'][stat] for stat in statList],width,color=colorDict['sbg'])
	# vglrRects=ax.bar(ind+7*width,[averageNormalizedStatDict['vglr'][stat] for stat in statList],width,color=colorDict['vglr'])
	# vg_haplo1kgRects=ax.bar(ind+8*width,[averageNormalizedStatDict['vg_haplo1kg'][stat] for stat in statList],width,color=colorDict['vg_haplo1kg'])

	# #add some text for labels, title and axes ticks
	# ax.set_title('Average, Normalized General Statistics')
	# ax.set_xticks(ind+width*4.5)
	# ax.set_xticklabels( statList )
	
	# box=ax.get_position()
	# ax.set_position([box.x0,box.y0,box.width*0.7,box.height])
	# ax.legend( (cactusRects[0],camelRects[0], k31Rects[0], k63Rects[0], prgRects[0],sbgRects[0],vglrRects[0],vg_haplo1kgRects[0]), ('Cactus','Camel','Debruijn k31','Debruijn k63','PRG','Seven Bridges','vglr','vg_haplo1kg'),loc='center left',bbox_to_anchor=(1, 0.5) )
	# plt.savefig('newStats.png')


	#Also output a tsv for each stat, where columns are regions and rows are algorithms
	for stat in statList:
		with open('tables/new/'+stat+'.tsv','w') as outFile:
			#Write header line
			outFile.write('\t'+'\t'.join([region.upper().replace('_',' ') for region in regionList])+'\n')
			for algo in sorted(statDict.keys(),key=lambda n:nameDict[n]):
				outFile.write(nameDict[algo]+'\t')
				for region in regionList:
					outFile.write(str(statDict[algo][region][stat])+'\t')
				outFile.write('\n')


if __name__=="__main__":
	main()