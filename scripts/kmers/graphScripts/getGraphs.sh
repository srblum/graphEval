#!/bin/sh
set -ex

#First argument should be url (with version number)
#Second argument should be output vg file like vg/new/BRCA1/cactus-brca1.vg


sg2vg $1 -u | vg view -Jv - > $2