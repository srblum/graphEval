#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from math import ceil
from collections import defaultdict
import matplotlib.pyplot as plt

def getData(inFile):
	dataDict=defaultdict(dict)
	with open(inFile) as inFile:
		for line in inFile:
			if not line.startswith('#'):
				line=line.strip()
				if line.startswith('$'):
					region=line[1:]
				else:
					line=line.split('\t')
					algo=line[0]
					if 'debruijn' in algo:
						if '31' in algo:
							algo='debruijn_k31'
						else:
							algo='debruijn_k63'
					else:
						algo=algo.split('-')[0]
					precision=line[2]
					recall=line[3]
					dataDict[region][algo]=[float(precision),float(recall)]
	return dataDict

def sortAlgoList(algoIndexPairs,nameDict):
	return sorted(algoIndexPairs,key=lambda n:nameDict[n[1]])


def main():
	inFile="data/completeThreshold1Results.txt"
	dataDict=getData(inFile)


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
		'trivial':'#b1b300',
		'snp1kg':'#fb9a99',
		'vglr':'#cab2d6',
		'haplo1kg30':'#00ff00',
		'haplo1kg50':'#0000ff',
		'shifted1kg':'#ff0000'
	}

	nameDict={
		'curoverse':'Curoverse',
		'debruijn_k31':'De Bruijn 31',
		'debruijn_k63':'De Bruijn 63',
		'cactus':'Cactus',
		'camel':'Camel',
		'sbg':'7BG',
		'refonly':'Primary',
		'trivial':'Unmerged',
		'vglr':'VGLR',
		'snp1kg':'1KG',
		'haplo1kg30':'1KG Haplo 30',
		'haplo1kg50':'1KG Haplo 50',
		'simons':'SGDP',
		'prg':'PRG',
		'shifted1kg':'Scrambled'
	}

	for region in ['BRCA1','SMA','LRC_KIR']:
		fig, ax = plt.subplots()
		algoList=[]
		scatterList=[]
		xMin=min(pair[0] for pair in dataDict[region].values())//0.01 * 0.01
		yMin=min(pair[1] for pair in dataDict[region].values())//0.01 * 0.01
		xMax=ceil(max(pair[0] for pair in dataDict[region].values())/0.01) * 0.01
		yMax=ceil(max(pair[1] for pair in dataDict[region].values())/0.01) * 0.01
		for algo in sorted(dataDict[region]):
			exec(algo+"=plt.scatter(dataDict[region][algo][0],dataDict[region][algo][1],c=[colorDict[algo]],alpha=0.8,s=[100],clip_on=False)")
			exec("algoList.append(algo)")
			exec("scatterList.append(eval(algo))")
		# trivialLine=plt.plot([0,1],[dataDict[region]['trivial'][1],dataDict[region]['trivial'][1]],color='#b1b300',ls='--')
		box=ax.get_position()
		ax.set_position([box.x0,box.y0,box.width,box.height])
		algoIndexPairs=zip(range(len(algoList)),algoList)
		sortedAlgoIndexPairs=sortAlgoList(algoIndexPairs,nameDict)
		sortedIndices=[pair[0] for pair in sortedAlgoIndexPairs]
		if region in ['BRCA1','SMA','LRC_KIR']:
			ax.set_position([box.x0,box.y0,box.width*0.7,box.height])
			plt.legend(tuple([scatterList[index] for index in sortedIndices]),tuple([nameDict[pair[1]] for pair in sortedAlgoIndexPairs]),scatterpoints=1,loc='upper left',fontsize=12,bbox_to_anchor=(1, 0.9))
		if region=='LRC_KIR': region='LRC KIR'
		ax.set_title(region,fontsize=22)
		ax.set_xlabel('Precision',fontsize=18)
		ax.set_ylabel('Recall',fontsize=18)

		#Get rid of the right and top borders
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')

		plt.ylim([yMin,1])
		plt.xlim([xMin,1])

		#Annotate all points in BRCA2 and MHC
		if region=='BRCA2':
			ax.annotate("Scrambled",(0.61,0.9835),color=colorDict['shifted1kg'],fontsize=14)
			ax.annotate("1KG",(0.83,0.996),color=colorDict['snp1kg'],fontsize=14)
			ax.annotate("Curoverse",(0.95,0.9965),color=colorDict['curoverse'],fontsize=14)
			ax.annotate("SGDP",(0.9,0.986),color=colorDict['simons'],fontsize=14)
			ax.annotate("Primary",(0.945,0.983),color=colorDict['refonly'],fontsize=14)
			ax.annotate("1KG Haplo 50",xy=(dataDict['BRCA2']['haplo1kg50'][0],dataDict['BRCA2']['haplo1kg50'][1]),xytext=(0.88,0.995),color=colorDict['haplo1kg50'],fontsize=14)
			ax.annotate("VGLR",xytext=(0.94,0.992),xy=(dataDict['BRCA2']['vglr'][0],dataDict['BRCA2']['vglr'][1]),color=colorDict['vglr'],arrowprops=dict(arrowstyle='-',color=colorDict['vglr']),fontsize=14)
			ax.annotate("PRG",xytext=(0.97,0.993),xy=(dataDict['BRCA2']['prg'][0],dataDict['BRCA2']['prg'][1]),color=colorDict['prg'],arrowprops=dict(arrowstyle='-',color=colorDict['prg']),fontsize=14)
			ax.annotate("De Bruijn 63",xytext=(0.88,0.990),xy=(dataDict['BRCA2']['debruijn_k63'][0],dataDict['BRCA2']['debruijn_k63'][1]),color=colorDict['debruijn_k63'],arrowprops=dict(arrowstyle='-',color=colorDict['debruijn_k63']),fontsize=14)
			ax.annotate("Unmerged",xytext=(0.90,0.989),xy=(dataDict['BRCA2']['trivial'][0],dataDict['BRCA2']['trivial'][1]),color=colorDict['trivial'],arrowprops=dict(arrowstyle='-',color=colorDict['trivial']),fontsize=14)
			ax.annotate("Cactus",xytext=(0.94,0.988),xy=(dataDict['BRCA2']['cactus'][0],dataDict['BRCA2']['cactus'][1]),color=colorDict['cactus'],arrowprops=dict(arrowstyle='-',color=colorDict['cactus']),fontsize=14)
			ax.annotate("Camel",xytext=(0.96,0.987),xy=(dataDict['BRCA2']['camel'][0],dataDict['BRCA2']['camel'][1]),color=colorDict['camel'],arrowprops=dict(arrowstyle='-',color=colorDict['camel']),fontsize=14)

		if region=='MHC':
			ax.annotate("Scrambled",(0.53,0.942),color=colorDict['shifted1kg'],fontsize=14)
			ax.annotate("1KG",(0.73,0.984),color=colorDict['snp1kg'],fontsize=14)
			ax.annotate("SGDP",(0.93,0.955),color=colorDict['simons'],fontsize=14)
			ax.annotate("Primary",(0.93,0.94),color=colorDict['refonly'],fontsize=14)
			ax.annotate("Camel",(0.56,0.984),color=colorDict['camel'],fontsize=14)
			ax.annotate("Cactus",(0.82,0.984),color=colorDict['cactus'],fontsize=14)

			ax.annotate("1KG Haplo 50",xy=(dataDict['MHC']['haplo1kg50'][0],dataDict['MHC']['haplo1kg50'][1]),xytext=(0.84,0.979),color=colorDict['haplo1kg50'],fontsize=14)
			
			ax.annotate("PRG",xytext=(0.97,0.99),xy=(dataDict['MHC']['prg'][0],dataDict['MHC']['prg'][1]),color=colorDict['prg'],arrowprops=dict(arrowstyle='-',color=colorDict['prg']),fontsize=14)
			ax.annotate("De Bruijn 63",xytext=(0.89,0.986),xy=(dataDict['MHC']['debruijn_k63'][0],dataDict['MHC']['debruijn_k63'][1]),color=colorDict['debruijn_k63'],arrowprops=dict(arrowstyle='-',color=colorDict['debruijn_k63']),fontsize=14)
			ax.annotate("Unmerged",xytext=(0.91,0.984),xy=(dataDict['MHC']['trivial'][0],dataDict['MHC']['trivial'][1]),color=colorDict['trivial'],arrowprops=dict(arrowstyle='-',color=colorDict['trivial']),fontsize=14)
			
		plt.savefig('newest_figures/'+region+'_complete_thresholded_weighted.png',format='png',dpi=600)

if __name__=="__main__":
	main()