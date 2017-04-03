This repo contains some of the code used to generate the data and visualizations present in the following paper on genome graphs:

http://biorxiv.org/content/early/2017/01/18/101378

Since the paper uses five genomic regions to evaluate graphs (BRCA1, BRCA2, SMA, MHC, and LRC_KIR), most of the data in this repo is divided into directories corresponding to those regions.  Graphs used during the analysis, as well as their urls on our temporary server, are listed in the file graph_servers.tsv in the base directory.

All code is located in the scripts directory.  

<strong>Please note that the directory structure of this repo was recently changed, and many of the filepaths hardcoded in the scripts need to be changed.</strong>

Code in the scripts directory is split into the following categories:

<h2>controlGraph</h2>
Generates a scrambled version of a graph, used as a control in our analyses.

<em>Requires vcf-sort, vg, and vg2sg</em>

<h2>cutDepth</h2>
Calls vg to perform a topological sort on a graph, and then cutDepth.py computes the average cut depth over the entire graph.  Note that vg must arbitrarily remove cycles to perform the sort.

<em>Requires vg</em>

<h2>kmers</h2>
Contains code to enumerate all 20mers in a set of fastq files--in our case, 2,691 low-coverage samples taken from the 1000 Genomes Project--as well as all 20mers in a set of graphs.  Then, for each graph, computes the portion of read-kmer instances present in the graph ("recall"), as well as the portion of graph-kmer instances present in the reads ("precision").

<em>Requires jellyfish and vg</em>

<h2>stats</h2>
Given the url of a ga4gh server graph, graphEval.py can perform a number of evaluations.  For more information, type:
```
graphEval.py -h
```
<h2>viz</h2>
Generates figures for the publication. 

<em>Requires python visualization library matplotlib.</em>