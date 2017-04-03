#!/usr/bin/env python
from __future__ import division,print_function
import sys


def checkForN(inFile):
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip()
			if line and len(line.split())==2:
				line=line.split()
				kmer=line[0]
				if 'N' in kmer:
					print('Found an n!')

def main():
	inFile=sys.argv[1]

	checkForN(inFile)

if __name__=="__main__":
	main()