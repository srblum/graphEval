#!/bin/sh
set -ex

#First argument should be input (.rev.fa file)
#Second argument should be output (.jf text file)

jellyfish count -m 20 -s 100 --text -o $2 $1