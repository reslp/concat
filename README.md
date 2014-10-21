concat.py
=========

python script to create concatenated FASTA files for phylogenetic analyses from single locus alignments


DESCRIPTION
===========

This script creates concatenated multilocus alignment files from single locus FASTA files for phylogenetic analyses.
It is under active development and may contain bugs. Please let me know when you use the script and tell me about your experiences.

How does it work?
================

`concat` takes multiple unaligned single locus FASTA files and a file containing the desired set of taxa as input.
It creates new single locus files reduced to the desired set of taxa, aligns the files and replaces gaps (-) at the beginning
and end of the alignments with question marks (?). Then it creates a concatenated alignment file on the basis of the aligned
single locus files and adds question marks for missing loci. All steps can also be executed individually.


REQUIREMENTS
============

- MacOS X or other Unix like operating system (Windows Version in the works)
- [python](www.python.org) 2.7.8+, which comes with most Unix like systems
- [mafft](http://mafft.cbrc.jp/alignment/software/) v7, for the alignment function


EXAMPLES
========

probably the simplest way to call concat.py is to provide a sequence ID file and a directory containing FASTA files of individual loci:

`python concat.py -t SeqIDFile.txt -d /path/to/sequences/ `

Getting help (displays all available command options):

`python concat.py -h`



COPYRIGTH AND LICENSE
=====================

Copyright (C) 2014 Philipp Resl

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program in the file LICENSE. If not, see http://www.gnu.org/licenses/.
