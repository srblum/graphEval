#!/usr/bin/env python
from __future__ import division, print_function
import argparse, string

def getRevComp(inputSeq):
	revSeq=inputSeq[::-1]
	transTable=string.maketrans('actgn','tgacn')
	revCompSeq=string.translate(revSeq,transTable)
	return revCompSeq

def parseArgs():
	parser = argparse.ArgumentParser(description="""Performs the specified evaluation on a specified graph server.
		Requires a url to be supplied from the user.""")
	type=parser.add_mutually_exclusive_group(required=True)
	parser.add_argument('--input',help="""The fasta file containing the sequences to be reverse complemented.""")
	parser.add_argument('--output',help="""The name of the output file containing both forward and reverse sequences.""")
	type.add_argument('--graph',action='store_true',default=False,help="""Specifies that the input file should be treated as a vg kmers output.""")
	type.add_argument('--fasta',action='store_true',default=False,help="""Specifies that the input file should be treated as a fasta file.""")
	args = parser.parse_args()
	return args


def getSequencesFromFasta(inFile):
	sequenceList=[]
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip().lower()
			if line and not line.startswith('>'):
				sequenceList.append(line)
				sequenceList.append(getRevComp(line))
	return sequenceList

def getSequencesFromVGKmers(inFile):
	sequenceList=[]
	with open(inFile) as inFile:
		for line in inFile:
			line=line.strip().lower().split()
			if line:
				kmer=line[0]
				sequenceList.append(kmer)
				sequenceList.append(getRevComp(kmer))
	return sequenceList

def writeSequences(outFile,sequenceList):
	with open(outFile,'w') as outFile:
		for sequence in sequenceList:
			outFile.write('>\n'+sequence+'\n')

def main():
	args=parseArgs()
	if args.graph:
		print(args)
		sequenceList=getSequencesFromVGKmers(args.input)
		writeSequences(args.output,sequenceList)
	elif args.fasta:
		sequenceList=getSequencesFromFasta(args.input)
		writeSequences(args.output,sequenceList)


if __name__=="__main__":
	main()
