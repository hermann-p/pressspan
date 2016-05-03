# Introduction

*pressspan* is a tool to visualise genomic data.  It has been
developed with data, that has been processed with
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) and
the *split* option, in mind. It can find reads mapped to different
chromosome strands, multiple chromosomes or which are circular.

The program uses massive multiprocessing to reassemble the genome and
identify such input reads.  During the assembly, all links containg a
specific fragment are merged into a tree-like graph.

Output can be filtered by a number of criteria. All elements matching
the given criteria are converted to [*.dot*
files](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29).

In addition, the Python script /select-subgraphs.py/ can be used to
pick results found on specific chromosomal positions and optionally
automatically create plots using your local /GraphViz/ packagage.

# Getting *pressspan*

## Binary version

The binary *presspan-xxx-standalone.jar* file is independent of system
and needs only a somewhat recent Java runtime environment (major
version 7 or later).

## Compile it

To build a snapshot version yourself, you will need
[Leiningen](http://leiningen.org/) and a proper Java SDK installed.

1) **get the source:** copy the "clone" link of the repository and execute
```
git clone (link)
```

2) **compile:**
```
cd pressspan
lein uberjar
```

3) When compiling is done, you will find a standalone version in
pressspan/target/uberjar/. Move it to a convenient place. You might
want to run
```
lein clean
```
when you are done to save disk space

# Usage

```
java -jar path/to/pressspan-xxx-standalone.jar -i (input file name) [options]
```

An input file is mandatory. It is advised to give at least one of
*-c*, *-m*, else reassembly is pointless.

Possible further options are:

## Basic options

**-h or --help** Show help and exit.

**-i or --in (path-to-infile)** Specify an input file. Currently only
   *unsorted .sam files created by segemehl* are supported and tested.

**-o or --out (path-to-outdir)** Specify a directory path for the
   output. If it doesn't exist, it gets created, if it does exist it
   does not become cleared. If none is given, *./pressspan_analysis*
   is used.

**-m or --multistrand** Detect reads containing fragments mapped to
   multiple chromosomes and write them as *.dot* files to
   *(out)/multis*

**-c or --circular** Detect reads containing fragments showing clues
   for circularity and write their *.dot* files to *(out)/circulars*

**-s or --bins** Adapt number of bins for frequency reports. Defaults
   to 10k bins per chromosome.

**-S or --stats** Write frequency report files to disk as .csv files.

## Filtering options

**-r or --range (chr:lower:upper)** Output only reads that contain at
   least one element the chromosome with ID-string (chr) between the
   positions (lower) and (upper). Does not (yet?) discriminate between
   plus and minus strand.

**-d or --drop (min-depth)** Don't output circular subgraphs in which
   any link's read depth is lower than (min-depth). Defaults to 1.

**-t or --trunc (min-depth)** When building multistrand subgraphs, do
   not include elements which are linked with a depth lower than
   (min-depth), truncate the graph there instead. Defaults to 1.

# Output

## Subgraphs

For each subgraphs matching your criteria, a *.dot* graph description
file is created in its corresponding subfolder. To visualise the
results you can use tools like [GraphViz](http://www.graphviz.org).

## Frequency statistics

The frequency of events found is reported to the tab-seperated files
*(out)/multistrand.csv* and *(out)/circulars.csv*. Within these files,
the first line contains the chromosome names, the second line contains
the size of each bin on the chromosome in bases, and all following
lines contain the number of times an event was detected in the
respective segment (that is, between (line-number - 2) * bin-size and
(line-nmuber - 1) * bin-size).

## Summary

In addition, a report file *(out)/pressspan.log* is written. It is a
tab-separated file with three columns:

1) **identifier string** which is also the name of the *.dot* file in
   its specific subdirectory

2) **list of fragment positions**, comma separated, in the format
   (strand-direction-prefix)(chromosome-id):(5'-position)-(3'-position)

3) **list of chromosomes** to which any reads have been matched,
   comma-separated
