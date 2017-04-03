#!/usr/bin/env python
from __future__ import division,print_function
import sys
import numpy as np
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
        'trivial':'#ffff99',
        'snp1000g':'#fb9a99',
        'vglr':'#cab2d6',
        'vg_haplo1kg':'#fdbf6f'
    }

    cactusPrecision=[0.7327,0.985,0.974]
    cactusRecall=[0.755,0.994,0.9736]
    camelPrecision=[0.67,0.977,0.931]
    camelRecall=[0.563,0.977,0.955]
    prgPrecision=[0.828,0.995,0.8213]
    prgRecall=[0.707,0.968,0.6295]
    k31Precision=[0.659,0.9409,0.8717]
    k31Recall=[0.6276,0.8688,0.6022]
    k63Precision=[0.6841,0.9854,0.9498]
    k63Recall=[0.6247,0.8444,0.6195]
    sbgPrecision=[0.986143338259,0.9901,0.738065463373]
    sbgRecall=[0.740126813604,0.97,0.659271450683]

    fig, (ax1,ax2,ax3) = plt.subplots(1,3, sharey=True)


    cactusSMA=ax1.scatter(cactusPrecision[0],cactusRecall[0],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
    camelSMA=ax1.scatter(camelPrecision[0],camelRecall[0],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
    prgSMA=ax1.scatter(prgPrecision[0],prgRecall[0],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
    sbgSMA=ax1.scatter(sbgPrecision[0],sbgRecall[0],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
    k31SMA=ax1.scatter(k31Precision[0],k31Recall[0],c=colorDict['debruijn_k31'],s=[100],alpha=0.8,marker='o')
    k63SMA=ax1.scatter(k63Precision[0],k63Recall[0],c=colorDict['debruijn_k63'],s=[100],alpha=0.8,marker='o')

    cactusMHC=ax3.scatter(cactusPrecision[1],cactusRecall[1],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
    camelMHC=ax3.scatter(camelPrecision[1],camelRecall[1],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
    prgMHC=ax3.scatter(prgPrecision[1],prgRecall[1],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
    sbgMHC=ax3.scatter(sbgPrecision[1],sbgRecall[1],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
    k31MHC=ax3.scatter(k31Precision[1],k31Recall[1],c=colorDict['debruijn_k31'],s=[100],alpha=0.8,marker='o')
    k63MHC=ax3.scatter(k63Precision[1],k63Recall[1],c=colorDict['debruijn_k63'],s=[100],alpha=0.8,marker='o')

    cactusLRC=ax2.scatter(cactusPrecision[2],cactusRecall[2],c=colorDict['cactus'],s=[100],alpha=0.8,marker='o')
    camelLRC=ax2.scatter(camelPrecision[2],camelRecall[2],c=colorDict['camel'],s=[100],alpha=0.8,marker='o')
    prgLRC=ax2.scatter(prgPrecision[2],prgRecall[2],c=colorDict['prg'],s=[100],alpha=0.8,marker='o')
    sbgLRC=ax2.scatter(sbgPrecision[2],sbgRecall[2],c=colorDict['sbg'],s=[100],alpha=0.8,marker='o')
    k31LRC=ax2.scatter(k31Precision[2],k31Recall[2],c=colorDict['debruijn_k31'],s=[100],alpha=0.8,marker='o')
    k63LRC=ax2.scatter(k63Precision[2],k63Recall[2],c=colorDict['debruijn_k63'],s=[100],alpha=0.8,marker='o')

    box1=ax1.get_position()
    ax1.set_position([box1.x0, box1.y0+box1.height*0.15,
                 box1.width, box1.height*0.8])

    box2=ax2.get_position()
    ax2.set_position([box2.x0, box2.y0+box2.height*0.15,
                 box2.width, box2.height*0.8])

    box3=ax3.get_position()
    ax3.set_position([box3.x0, box3.y0+box3.height*0.15,
                 box3.width, box3.height*0.8])


    


    # ax.legend((sbgSMA,
    #     cactusSMA,camelSMA,
    #     k31SMA,k63SMA,
    #     prgSMA),
    # ('7BG','Cactus','Camel', 'De Bruijn 31','De Bruijn 63','PRG',),
    #     scatterpoints=1,loc='center left',bbox_to_anchor=(1,0.5))

    plt.figlegend((sbgMHC,k31MHC,cactusMHC,k63MHC,camelMHC,prgMHC),
    ('7BG','De Bruijn 31','Cactus','De Bruijn 63','Camel','PRG',),
        scatterpoints=1,loc='center left',ncol=3,fontsize=10,bbox_to_anchor=(0.25,0.08))

    # ax.legend((sbgLRC,
    #     cactusLRC,camelLRC,
    #     k31LRC,k63LRC,
    #     prgLRC),
    # ('7BG','Cactus','Camel', 'De Bruijn 31','De Bruijn 63','PRG'),
    #     scatterpoints=1,loc='center left',bbox_to_anchor=(1,0.5))




    # ax.legend((sbgSMA,sbgMHC,sbgLRC,
    #     cactusSMA,cactusMHC,cactusLRC,camelSMA,camelMHC,camelLRC,
    #     k31SMA,k31MHC,k31LRC,k63SMA,k63MHC,k63LRC,
    #     prgSMA,prgMHC,prgLRC),
    # ('7BG SMA','7BG MHC','7BG LRC KIR','Cactus SMA','Cactus MHC','Cactus LRC KIR',
    #     'Camel SMA','Camel MHC','Camel LRC KIR', 'De Bruijn 31 SMA', 'De Bruijn 31 MHC','De Bruijn 31 LRC KIR',
    #     'De Bruijn 63 SMA', 'De Bruijn 63 MHC','De Bruijn 63 LRC KIR','PRG SMA','PRG MHC', 'PRG LRC KIR',),
    #     scatterpoints=1,ncol=6,loc='center left',fontsize=7,bbox_to_anchor=(-0.15, -0.3))
    plt.setp( (ax1,ax2,ax3) ,xticks=[0.6,0.7,0.8,0.9,1.0])

    ax1.set_ylim([0.5,1])
    ax1.set_xlim([0.6,1])
    ax2.set_ylim([0.5,1])
    ax2.set_xlim([0.6,1])
    ax3.set_ylim([0.5,1])
    ax3.set_xlim([0.6,1])
    fig.suptitle('GRC Alignment Concordance')
    ax1.set_title('SMA')
    ax2.set_title('LRC')
    ax3.set_title('MHC')
    # ax1.set_xlabel('Precision')
    ax2.set_xlabel('Precision')
    # ax3.set_xlabel('Precision')
    ax1.set_ylabel('Recall')
    plt.savefig('GRC_Alignment_Subplots.png',format='png',dpi=600)
if __name__=="__main__":
    main()