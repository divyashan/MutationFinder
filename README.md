# MutationFinder
An open source library built for identifying mentions of mutations within bodies of text.

A regular expression driven library for extracting mutation mentions (in both protein and chromosome formats) from bodies of text. More details can be found [here](http://mutationfinder.sourceforge.net/).

### Running MutationFinder

If you have a file formatted as described in (Caporaso et al., 2007) and
the input file format section below, you can apply MutationFinder with 
the following steps:

> cd MutationFinder
> ./mutation_finder.py /path/to/your/input/file

A new file will be created in the current working directory called 
 input_filename.mf 

A non-default output directory can be specified with the -o flag. Run:

> ./mutation_finder.py -h

for more information on parameters which can be passed to mutation_finder.py

### INPUT FILE FORMAT

The input files to be processed by MutationFinder should contain one 'document'
per line. Each line should be tab-delimited and contain two fields: a document 
identifier and the document text. See the devo_set.txt and test_set.txt files in
the MutationFinder/corpora directory for examples.




MutationFinder: A high-performance system for extracting point mutation mentions from text
J. Gregory Caporaso, William A. Baumgartner Jr., David A. Randolph, K. Bretonnel Cohen, and Lawrence Hunter; Bioinformatics, 2007 23(14):1862-1865; doi:10.1093/bioinformatics/btm235;
