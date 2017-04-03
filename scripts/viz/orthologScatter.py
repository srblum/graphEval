#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def getData(inFile):
	def getZeroList():
		return [[0,0],[0,0],[0,0]]
	orthologDict={
		'lrc_kir':21199228,
		'sma':1020503,
		'mhc':10811041
	}
	baseDict={
		'lrc_kir':658780,
		'sma':748036,
		'mhc':2022557
	}
	altDict={
		'lrc_kir':35,
		'sma':2,
		'mhc':8
	}
	regionList=[]
	dataDict=defaultdict(getZeroList)
	regionCount=-1
	algo=''
	source=''
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line.startswith('$'):
				region=line[1:]
				regionList.append(region)
				regionCount+=1
			elif line.startswith('@'):
				line=line.split()
				algo=line[0][1:]
				try:
					source=line[1]
				except:
					print(line)
					sys.exit('Test complete.')
			else:
				line=line.split('\t')
				allele=line[0]
				try:
					orthologs=line[1]
					paralogs=line[2]
				except:
					print(line)
				if allele=='total':
					if algo!='trivial':
						dataDict[algo][regionCount]=[int(orthologs)/orthologDict[region],int(paralogs)/baseDict[region]/altDict[region]]
	return dataDict, regionList



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
		'trivial':'#ffff99',
		'snp1000g':'#fb9a99',
		'vglr':'#cab2d6',
		'vg_haplo1kg':'#fdbf6f'
	}

	inFile="noN_gene_results.txt"
	dataDict, regionList=getData(inFile)

	fig, (ax1,ax2,ax3) = plt.subplots(1,3, sharey=True)


	cactusSMA=ax1.scatter(dataDict['cactus'][0][0],dataDict['cactus'][0][1],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
	camelSMA=ax1.scatter(dataDict['camel'][0][0],dataDict['camel'][0][1],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
	prgSMA=ax1.scatter(dataDict['prg'][0][0],dataDict['prg'][0][1],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
	sbgSMA=ax1.scatter(dataDict['sbg'][0][0],dataDict['sbg'][0][1],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
	debruijn_k31SMA=ax1.scatter(dataDict['debruijn_31k'][0][0],dataDict['debruijn_31k'][0][1],c=colorDict['debruijn_31k'],s=[100],alpha=0.8,marker='o')
	debruijn_k63SMA=ax1.scatter(dataDict['debruijn_63k'][0][0],dataDict['debruijn_63k'][0][1],c=colorDict['debruijn_63k'],s=[100],alpha=0.8,marker='o')
	
	cactusMHC=ax3.scatter(dataDict['cactus'][1][0],dataDict['cactus'][1][1],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
	camelMHC=ax3.scatter(dataDict['camel'][1][0],dataDict['camel'][1][1],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
	prgMHC=ax3.scatter(dataDict['prg'][1][0],dataDict['prg'][1][1],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
	sbgMHC=ax3.scatter(dataDict['sbg'][1][0],dataDict['sbg'][1][1],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
	debruijn_k31MHC=ax3.scatter(dataDict['debruijn_31k'][1][0],dataDict['debruijn_31k'][1][1],c=colorDict['debruijn_31k'],s=[100],alpha=0.8,marker='o')
	debruijn_k63MHC=ax3.scatter(dataDict['debruijn_63k'][1][0],dataDict['debruijn_63k'][1][1],c=colorDict['debruijn_63k'],s=[100],alpha=0.8,marker='o')
	
	cactusLRC=ax2.scatter(dataDict['cactus'][2][0],dataDict['cactus'][2][1],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
	camelLRC=ax2.scatter(dataDict['camel'][2][0],dataDict['camel'][2][1],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
	prgLRC=ax2.scatter(dataDict['prg'][2][0],dataDict['prg'][2][1],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
	sbgLRC=ax2.scatter(dataDict['sbg'][2][0],dataDict['sbg'][2][1],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
	debruijn_k31LRC=ax2.scatter(dataDict['debruijn_31k'][2][0],dataDict['debruijn_31k'][2][1],c=colorDict['debruijn_31k'],s=[100],alpha=0.8,marker='o')
	debruijn_k63LRC=ax2.scatter(dataDict['debruijn_63k'][2][0],dataDict['debruijn_63k'][2][1],c=colorDict['debruijn_63k'],s=[100],alpha=0.8,marker='o')

	box1=ax1.get_position()
	ax1.set_position([box1.x0, box1.y0+box1.height*0.15,box1.width, box1.height*0.8])

	box2=ax2.get_position()
	ax2.set_position([box2.x0, box2.y0+box2.height*0.15,box2.width, box2.height*0.8])

	box3=ax3.get_position()
	ax3.set_position([box3.x0, box3.y0+box3.height*0.15,box3.width, box3.height*0.8])

	plt.figlegend((sbgMHC,debruijn_k31MHC,cactusMHC,debruijn_k63MHC,camelMHC,prgMHC),('7BG','De Bruijn 31','Cactus','De Bruijn 63','Camel','PRG',),
		scatterpoints=1,loc='center left',ncol=3,fontsize=10,bbox_to_anchor=(0.25,0.08))
	# plt.legend((sbgSMA,sbgMHC,sbgLRC,cactusSMA,cactusMHC,cactusLRC,camelSMA,camelMHC,camelLRC,\
	# 	debruijn_k31SMA,debruijn_k31MHC,debruijn_k31LRC,debruijn_k63SMA,debruijn_k63MHC,debruijn_k63LRC,\
	# 	prgSMA,prgMHC,prgLRC),\
	# 	('7BG SMA','7BG MHC','7BG LRC KIR','Cactus SMA','Cactus MHC','Cactus LRC KIR','Camel SMA','Camel MHC','Camel LRC KIR',\
	# 	'De Bruijn 31 SMA', 'De Bruijn 31 MHC','De Bruijn 31 LRC KIR','De Bruijn 63 SMA','De Bruijn 63 MHC','De Bruijn 63 LRC KIR',\
	# 	'PRG SMA','PRG MHC','PRG LRC KIR'),\
	# 	scatterpoints=1,ncol=6,loc='upper left',fontsize=7,bbox_to_anchor=(-0.15, -0.15))
	ax1.set_ylim([0,.03])
	ax1.set_xlim([0,1])
	ax2.set_ylim([0,.03])
	ax2.set_xlim([0,1])
	ax3.set_ylim([0,.03])
	ax3.set_xlim([0,1])

	plt.setp( (ax1,ax2,ax3) ,xticks=[0.2,0.4,0.6,0.8,1.0])

	ax1.set_title('SMA')
	ax2.set_title('LRC KIR')
	ax3.set_title('MHC')

	fig.suptitle('Ortholog vs. Paralog Mappings')
	ax2.set_xlabel('Ortholog Mappings per Possible Ortholog Mapping')
	ax1.set_ylabel('Paralog Mappings per Genic Base per Alt Path')
	plt.savefig('noN_Ortholog_Paralog_Subplots.png',format='png',dpi=600)

if __name__=="__main__":
	main()