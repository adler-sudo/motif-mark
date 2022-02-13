# Motif Mark Plan

Generate plan of attack for creating visual of motif locations within a series of sequences. 

## List of classes

* Read
    * takes each sequence from the fasta
    * some kind of unique identifier
    * keep track of which motifs match where
    * location of exons
    * location of introns

* Motif
    * contains motif sequence
    * functionality for ambiguous nucleotides

* Gene - not clear on read vs. gene?

## List of functions

* pycairo drawing function
* argparse
    * -f: fasta file
    * -m: motifs file
    * -o: output file (same prefix as fasta with .png)
* read through fasta
* read through motif file
* function to match motifs within reads/genes
