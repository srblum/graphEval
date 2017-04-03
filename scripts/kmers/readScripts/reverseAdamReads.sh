#!/bin/sh
set -e

fastx_reverse_complement -i $1 -o $2
cat $1 >> $2