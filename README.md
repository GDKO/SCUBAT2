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
Requiries a BLAST XML file i.e
```blastn -query transcripts.fa -db contigs.fa -evalue 1e-25 -outfmt 5 -out blast.xml```

