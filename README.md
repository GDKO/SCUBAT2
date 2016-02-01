SCUBAT2
===============
###Overview

SCUBAT2 (Scaffolding Contigs Using BLAST And Transcripts v2) uses transcriptome or proteome information to scaffold the genome. It was inspired by the original [SCUBAT](https://github.com/elswob/SCUBAT/) algorithm by Ben Elsworth.

Requirements
------------
###Python Libraries 
[Biopython](http://biopython.org/wiki/Main_Page) - to parse BLAST XML file

[Numpy](http://www.numpy.org/) - to calculate some statistics

Details
------------
Requires a BLAST XML file

```blastn -query transcripts.fa -db contigs.fa -evalue 1e-25 -outfmt 5 -out blast.xml```

For the same species the default settings for identity cutoff should be okay

The user must specify the max allowed intron size (i.e for nematode species ~ 20000 bp). Alternatively the user can run the program with ```--intron_size_run``` that creates the file intron_size which has the intron sizes calculated by the mapped transcripts