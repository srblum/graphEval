#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

colorDict={
	'curoverse':'#a6cee3',
	'debruijn-k31':'#e31a1c',
	'debruijn-k63':'#ff7f00',
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
	'debruijn-k31':'De Bruijn 31',
	'debruijn-k63':'De Bruijn 63',
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


def getData(inFile):
	dataDict={}
	with open(inFile) as inFile:
		for line in inFile:
			if not line.startswith('#'):
				line=line.strip().split('\t')
				algo=line[0].split('_')[0]
				total=float(line[1])
				nonRef=float(line[2])
				if algo not in ['trivial','sbg','camel']:
					dataDict[algo]=[total,nonRef]
	return dataDict


dataDir='calls/'
figDir='medianFigs/'
regionList=['brca2','mhc']
interval=10
for region in regionList:
	inFile=dataDir+region+'.tsv'
	dataDict=getData(inFile)
	fig, ax = plt.subplots(figsize=(16,6))
	algoList=[]
	tickList=[]
	# colorList=[]
	for index,pair in enumerate(sorted(dataDict.items(),reverse=True,key=lambda n:n[1][1])):
		algo, valuePair=pair
		algoList.append(nameDict[algo])
		tickList.append(index*interval+4)
		# colorList.append(colorDict[algo])
		if region=='mhc':
			total,nonRef=[value/1000 for value in valuePair]
		else:
			total,nonRef=valuePair
		blackBar=plt.bar(4+index*interval,total,3,color='black') #plot the total calls
		redBar=plt.bar(7+index*interval,nonRef,3,color='red') #plot the nonRef calls

	#Resize to fit labels and legend
	box=ax.get_position()
	ax.set_position([box.x0,box.y0+box.height*0.18,box.width*0.9,box.height*0.82])

	#Get rid of the right and top borders
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('left')

	#Make Legend
	ax.legend( (blackBar[0],redBar[0]),('# Calls','# Non-ref Calls'),loc='center left',bbox_to_anchor=(1, 0.5) )

	#Set title
	ax.set_title(region.upper(),fontsize=28)

	#Label y axis
	if region=='mhc':
		ax.set_ylabel('Calls (Thousands)',fontsize=24)
	else:
		ax.set_ylabel('Calls',fontsize=24)

	#Label x axis
	ax.set_xticks(tickList)
	ax.set_xticklabels(algoList,rotation=45,fontsize=18)
	# for label,color in zip(ax.get_xticklabels(),colorList):
	# 	label.set_color(color)

	# save
	plt.savefig(figDir+region+'_calls_noTrivial_no7BG_noCamel_smallTitle.pdf',format='pdf')
