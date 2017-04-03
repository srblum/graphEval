#!/usr/bin/env python
from __future__ import division, print_function
import argparse, sys
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy import stats

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

# Linear Regression
def polyfit(x, y):
    coeffs = np.polyfit(x, y, 1)

    # r-squared
    p = np.poly1d(coeffs)

    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    rSquared= round(ssreg / sstot,2)
    return rSquared

def getData(inFile):
	dataDict={}
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip().split('\t')
			algo=line[0].split('_')[0]
			try:
				vcf=float(line[3])
				kmer=float(line[6])
				corg=float(line[7])
				dataDict[algo]=[kmer,corg,vcf]
			except:
				pass
	return dataDict


dataDir='glenn/'
figDir='../../figs/variants/'
regionList=['brca1','brca2','mhc']



reverseNameDict={}

for x,y in nameDict.items():
	reverseNameDict[y]=x

pDict=defaultdict(list)

dataDict=defaultdict(list)
for region in regionList:
	inFile=dataDir+"fixed"+region+'.txt'
	regionDict=getData(inFile)
	for algo,list_ in regionDict.items():
		dataDict[algo].append(list_)



statList=["K-mer F1","1-Corg"]
fig, axarr = plt.subplots(1,2)
for iStat,stat in enumerate(statList):
	
	# axarr[0,iStat].set_title(stat)
	xList=[]
	yList=[]
	for algo in dataDict:
		for point in dataDict[algo]:
			xList.append(point[2])
			yList.append(point[iStat])
			axarr[iStat].scatter(point[2],point[iStat],color=colorDict[reverseNameDict[algo]])
	# Show best fit line
	m,b=np.polyfit(xList,yList,1)
	axarr[iStat].plot(xList,[m*x+b for x in xList],'-',color='black')

	# Show rSquared
	rSquared=polyfit(xList,yList)
	axarr[iStat].annotate("$R^2={}$".format(rSquared),xy=(0.02,0.02),xycoords='axes fraction')

	#Set stat labels (y axes)
	axarr[iStat].set_xlabel(stat)

	#Set fScore label (single x axis label)
	axarr[0].set_ylabel('VCF F1')

	#Get rid of ticks
	[ax.xaxis.set_ticks_position('none') for ax in axarr.flatten()]
	[ax.yaxis.set_ticks_position('none') for ax in axarr.flatten()]

	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in axarr.flatten()], visible=False)
	plt.setp([a.get_yticklabels() for a in axarr.flatten()], visible=False)

#Save
plt.savefig(figDir+"vcfF1vsStats.png",format='png',dpi=400)




	# # Compute spearmanr and convert to one-tailed
	# r,twoTailP=stats.spearmanr(xList,yList)
	# print(stat, r, twoTailP)

	# Compute pearsonr and convert to one-tailed
	# r,twoTailP=stats.pearsonr(xList,yList)
	# print(stat, r, twoTailP)




	# sxList=sorted(enumerate(xList),key=lambda n:n[1][0])
	# syList=sorted(enumerate(yList),key=lambda n:n[1][0])
	# sDict=defaultdict(dict)
	# for i,l in enumerate([sxList,syList]):
	# 	for rank,pair in enumerate(l):
	# 		id_,pair=pair
	# 		score,algo=pair
	# 		sDict[id_][i]=[rank,algo]
	# fig, ax = plt.subplots()
	# for id_ in sDict:
	# 	ax.scatter(sDict[id_][0][0],sDict[id_][1][0],color=colorDict[reverseNameDict[sDict[id_][0][1]]])
	# ax.set_title("VCF F1 rank vs. "+stat+" rank")
	# ax.set_xlabel('VCF F1 rank')
	# ax.set_ylabel(stat+" rank")
	# plt.savefig(figDir+"vcfF1vs"+stat+"Ranks.png""",format='png')
	# if stat in ['Compression','Base Degree']:
	# 	if r<=0:
	# 		print(stat,r,oneTailP)
	# 		oneTailP=1-oneTailP
	# else:
	# 	if r>=0:
	# 		print(stat,r,oneTailP)
	# 		oneTailP=1-oneTailP
	# pDict[stat].append(oneTailP)

		#Show spearman r
# 		axarr[iStat,iRegion].annotate("$r={0:.3f}$".format(r),xy=(0.02,0.02),xycoords='axes fraction')

# for stat,list_ in pDict.items():
# 	combinedP=stats.combine_pvalues(list_)
# 	print("Stat: {}\nP value: {}".format(stat,combinedP))


