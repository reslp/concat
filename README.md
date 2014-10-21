concat.py
=========

python script to create concatenated FASTA files for phylogenetic analyses from single locus alignments


DESCRIPTION
===========

This script creates concatenated multilocus alignment files from single locus FASTA files for phylogenetic analyses.
It is under active development and may contain bugs. Please let me know when you use the script and tell me about your experiences.

How does it work?
================

concat takes multiple unaligned single locus files and a file containing the desired set of taxa as input.
It creates new single locus files reduced to the desired set of taxa, aligns the files and replaces gaps (-) at the beginning
and end of the alignments with question marks (?). Then it creates a concatenated alignment file an the basis of the aligned
single locus files and adds question marks for missing loci.

