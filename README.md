# Introduction

*pressspan* is a tool to visualise genomic data.
It has been developed with data, that has been processed with
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/)
and the *split* option, in mind. It can find reads mapped to different chromosome strands, multiple chromosomes or which are circular.

The program uses massive multiprocessing to reassemble the genome and identify such input reads.
During the assembly, all links containg a specific fragment are merged into a tree-like graph.

Output can be filtered by a number of criteria. All elements matching the given criteria are converted to [*.dot* files](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29).

# Getting *pressspan*

## Binary version

The binary *presspan-xxx-standalone.jar* file is independent of system and needs only a somewhat recent Java runtime environment (tested on Java major versions 7 and 8).

## Compile it

To build a snapshot version yourself, you will need [Leiningen](http://leiningen.org/) and a proper Java SDK installed.

1) **get the source:** copy the "clone" link of the repository and execute
```
git clone (link)
```

2) **compile:**
```
cd pressspan
lein uberjar
```

3) When compiling is done, you will find a standalone version in pressspan/target/uberjar/. Move it to a convenient place. You might want to run
```
lein clean
```
when you are done to save disk space

# Usage

```
java -jar path/to/pressspan-xxx-standalone.jar -i (input file name) [options]
```

An input file needs to be given. Possible options are:

## Basic options

**-h or --help** Show help and exit.

**-i or --in (path-to-infile)** Specify an input file. Currently only *unsorted .sam files created by segemehl* are supported and tested.

**-o or --out (path-to-outdir)** Specify a directory path for the output. If it doesn't exist, it gets created, if it does exist it does not become cleared. If none is given, *./pressspan_analysis* is used.

**-m or --multistrand** Detect reads containing fragments mapped to multiple chromosomes and write them as *.dot* files to *(out)/multis*

**-c or --circular** Detect reads containing fragments showing clues for circularity and write their *.dot* files to *(out)/circulars*

## Filtering options

**-r or --report (path-to-file)** Give filename for the report file. OBSOLETE. Defaults to *(out)/pressspan.log*.

**-R or --range (chr:lower:upper)** Output only reads that contain at least one element the chromosome with ID-string (chr) between the positions (lower) and (upper). Does not (yet?) discriminate between plus and minus strand.

**-d or --drop (min-depth)** Don't output subgraphs in which any link's read depth is lower than (min-depth). Defaults to 1.

**-t or --trunc (min-depth)** When building subgraphs, do not include elements which are linked with a depth lower than (min-depth), truncate the graph there instead. Defaults to 1.

**-p or --position (chr:p5-pos)** Find all possible links to a given fragment. NOT IMPLEMENTED YET.

# Output

For each subgraphs matching your criteria, a *.dot* graph description file is created in its corresponding subfolder. To visualise the results you can use tools like [GraphViz](http://www.graphviz.org).

In addition, a report file *(out)/pressspan.log* is written. It is a tab-separated file with three columns:

1) **identifier string** which is also the name of the *.dot* file in its specific subdirectory

2) **list of fragment positions**, comma separated, in the format (strand-direction-prefix)(chromosome-id):(5'-position)-(3'-position)

3) **list of chromosomes** to which any reads have been matched, comma-separated
