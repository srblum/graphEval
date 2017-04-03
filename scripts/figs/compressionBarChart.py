#!/usr/bin/env python
# a bar plot with errorbars
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt



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
		'trivial':'#b1b300',
		'snp1kg':'#fb9a99',
		'vglr':'#cab2d6',
		'haplo1kg30':'#00ff00',
		'haplo1kg50':'#0000ff',
		'shifted1kg':'#ff0000'
	}
	def getData(inFile,regionList):
		def getZeroList():
			return [0,0,0,0,0]
		dataDict=defaultdict(getZeroList)
		algo=''
		with open(inFile) as inFile:
			for line in inFile:
				line=line.strip()
				if line.startswith('$'):
					region=line[1:]
				else:
					if region!='cenx':
						line=line.split('\t')
						algo=line[0]
						positions=int(line[2])
						nonNPositions=line[3]
						baseJoins=float(line[8])
						# dataDict[algo][regionList.index(region)]=nonNRefPositionDict[region]/int(nonNPositions)
						dataDict[algo][regionList.index(region)]=round(baseJoins/positions,4)

		return dataDict

	#Dict of non-N primary reference length for each region
	#Note that there are no N's in the primary reference anyways...
	nonNRefPositionDict={
		'sma':2397625,
		'lrc_kir':1058685,
		'mhc':4970458,
		'brca1':81189,
		'brca2':84989
	}
	inFile="base_degree_stats_results.txt"
	regionList=['brca1','brca2','lrc_kir','mhc','sma']
	dataDict=getData(inFile,regionList)
	N=5
	ind = np.arange(N)*2+0.2  # the x locations for the groups
	width = 0.11       # the width of the bars

	fig, ax = plt.subplots()
	print(dataDict)
	kgRects=ax.bar(ind+width,dataDict['snp1kg'],width,color=colorDict['snp1kg'])
	haplo30Rects=ax.bar(ind+2*width,dataDict['haplo1kg30'],width,color=colorDict['haplo1kg30'])
	haplo50Rects=ax.bar(ind+3*width,dataDict['haplo1kg50'],width,color=colorDict['haplo1kg50'])
	sbgRects=ax.bar(ind+4*width,dataDict['sbg'],width,color=colorDict['sbg'])
	cactusRects=ax.bar(ind+5*width,dataDict['cactus'],width,color=colorDict['cactus'])
	camelRects=ax.bar(ind+6*width,dataDict['camel'],width,color=colorDict['camel'])
	controlRects=ax.bar(ind+7*width,dataDict['shifted1kg'],width,color=colorDict['shifted1kg'])
	curoRects=ax.bar(ind+8*width,dataDict['curoverse'],width,color=colorDict['curoverse'])
	k31Rects=ax.bar(ind+9*width,dataDict['debruijn_31k'],width,color=colorDict['debruijn_k31'])
	k63Rects=ax.bar(ind+10*width,dataDict['debruijn_63k'],width,color=colorDict['debruijn_k63'])
	prgRects=ax.bar(ind+11*width,dataDict['prg'],width,color=colorDict['prg'])
	simonRects=ax.bar(ind+12*width,dataDict['simons'],width,color=colorDict['simons'])
	vglrRects=ax.bar(ind+13*width,dataDict['vglr'],width,color=colorDict['vglr'])

	#add some text for labels, title and axes ticks
	ax.set_ylabel('Base Degree')
	# ax.set_title('Compression (Primary Reference Positions/Non-N Graph Positions)')
	ax.set_title('Base Degree')
	ax.set_xticks(ind+width*4)
	ax.set_xticklabels(['BRCA1','BRCA2','LRC KIR','MHC','SMA'])
	
	box=ax.get_position()
	ax.set_position([box.x0,box.y0,box.width*0.75,box.height])
	ax.legend(
		(kgRects[0],haplo30Rects[0],haplo50Rects[0],sbgRects[0],cactusRects[0],camelRects[0],controlRects[0],
			curoRects[0],k31Rects[0], k63Rects[0], prgRects[0],simonRects[0],vglrRects[0]), 
		('1KG','1KG Haplo 30','1KG Haplo 50','7BG','Cactus','Camel','Control',
			'Curoverse','De Bruijn 31','De Bruijn 63','PRG','SGDP','VGLR'),loc='center left',bbox_to_anchor=(1, 0.5) )
	plt.savefig('baseDegree.png',format='png',dpi=600)



if __name__=="__main__":
	main()