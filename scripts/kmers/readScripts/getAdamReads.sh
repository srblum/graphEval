#!/bin/sh
set -e

#First argument should be the path to a fastq file.
#Second argument should be the path to a fasta file.

fastq_to_fasta -Q 33 -i $1 -o $2