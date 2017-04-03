#!/bin/sh
set -e

#First argument should be vg file like vg/new/BRCA1/curoverse.vg
#Second argument should be output file like kmers/new/BRCA1/curoverse_edgeLimit7.kmer
vg kmers -k 20 -e 7 -t 8 -p $1 > $2
