This repo contains some of the code used to generate the data and visualizations present in the following paper on genome graphs:

http://biorxiv.org/content/early/2017/01/18/101378

Code in the scripts directory is split into the following categories:

<h2>controlGraph</h2>
Generates a scrambled version of a graph, used as a control in our analyses.
<br><br>
<em>Requires vcf-sort, vg, and vg2sg</em>

<h2>cutDepth</h2>
Calls vg to perform a topological sort on a graph, and then cutDepth.py computes the average cut depth over the entire graph.  Note that vg must arbitrarily remove cycles to perform the sort.
<br><br>
<em>Requires vg</em>

<h2>figs</h2>
Generates figures for the paper. 
<br><br>
<em>Requires the python visualization library matplotlib.</em>

<h2>graphViz</h2>
Creates a dynamic browser visualization using data from a user-specified ga4gh graph server url.  (Demo currently unavailable because server is down.)

<h2>kmers</h2>
Contains code to enumerate all kmers in a set of fastq files--in our case, 2,691 low-coverage samples taken from the 1000 Genomes Project--as well as all kmers in a set of graphs.  Then, for each graph, computes the portion of read-kmer instances present in the graph ("recall"), as well as the portion of graph-kmer instances present in the reads ("precision").
<br><br>
<em>Requires jellyfish and vg</em>

<h2>stats</h2>
Given the url of a ga4gh server graph, graphEval.py can perform a number of evaluations.  See the documentation in graphEval.py for details.

