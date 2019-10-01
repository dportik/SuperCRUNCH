![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Overview

**SuperCRUNCH** is a python toolkit for creating and working with phylogenetic datasets. SuperCRUNCH can be run using any set of sequence data, as long as sequences are in fasta format with standard naming conventions (described [here](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Sequences)). 

SuperCRUNCH can be used to process sequences downloaded directly from GenBank/NCBI, local sequence data (e.g. sequences not downloaded from GenBank, such as unpublished data), or a combination of both. The sequence data are first parsed into gene-specific fasta files using targeted searches guided by lists of taxon and locus names. For each resulting gene, sequences can be filtered with similarity searches using automated methods or based on user-supplied reference sequences. SuperCRUNCH offers the option to select a best representative sequence for each taxon, or to retain all filtered sequences for each taxon. These options allow the user to generate species-level supermatrix datasets (one sequence per species per locus) or population-level datasets (multiple sequences per species per locus). In addition, SuperCRUNCH can identify voucher codes present in sequence records and link samples in phylogeographic datasets through correct labeling. SuperCRUNCH offers important pre-alignment steps (adjust sequence directions, adjust reading frames), several options for sequence alignment (Clustal-O, MAFFT, Muscle, MACSE), and multiple options for alignment trimming. Finally, SuperCRUNCH can be used for rapid file format conversion and concatenation. 

SuperCRUNCH is scalable and can be used to assemble a variety of datasets, ranging from small population-level datasets (one taxon, one gene) to large phylogenomic datasets with thousands of loci (such as UCEs or other sequence capture datasets). SuperCRUNCH was intended to be transparent, objective and repeatable, and provides meaningful output at every step to help guide user decisions. In addition, it is modular in design and various components of SuperCRUNCH can be easily incorporated into other custom bioinformatics workflows.

A full overview of SuperCRUNCH is described in the following pre-print article:

+ Portik, D.M., and J.J. Wiens. (2019) SuperCRUNCH: A toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. BioRxiv, https://doi.org/10.1101/538728.


## Version

The current release of **SuperCRUNCH** is [**v1.2**](https://github.com/dportik/SuperCRUNCH/releases). 

### Major changes in v1.2:
  - **All modules are now compatible with Python 2.7 and Python 3.7.**
  - **SQL now implemented for creating local databases and performing searches.** Improved modules include `Parse_Loci.py` (~30x speedup), `Filter_Seqs_and_Species.py` (~20x speedup), `Taxon_Assessment.py` (~3x speedup), and `Rename_Merge.py` (~3x speedup).
  - **Ability to identify and use voucher information to create 'vouchered' datasets.** This feature can link population-level samples across loci, allowing rapid reconstruction of published phylogeographic datasets.
  - Added output directory specification to all modules.
  - Two trimming modules now included: `Trim_Alignments_Trimal.py` and `Trim_Alignments_Custom.py`. The `Trim_Alignments_Custom.py` module allows finding start and stop block positions, and row-wise (internal) sliding window trimming based on divergence.
  - Added new module `Filter_Fasta_by_Min_Seqs.py` to filter fasta files using a minimum number of sequences.
  - Added multithreading for BLAST searches and new `--bp_bridge` flag for coordinate merging in `Cluster_Blast_Extract.py` and `Reference_Blast_Extract.py`.
  - Added `--multisearch` option to `Reference_Blast_Extract.py` and `Contamination_Filter.py` to automate sequential runs within a directory.
  - Complete code re-write for several modules, including `Align.py`, `Cluster_Blast_Extract.py`, `Filter_Seqs_and_Species.py`, `Parse_Loci.py`, `Taxon_Assessment.py`.
  - Re-ordered tasks in `Cluster_Blast_Extract.py` to allow completion of all steps for one fasta file before moving to next fasta file in sequence.
  - Module `Relabel_Fasta.py` was renamed to `Fasta_Relabel_Seqs.py`.

For complete version history please see the [change log file](https://github.com/dportik/SuperCRUNCH/tree/master/CHANGELOG.md).


## Installation

SuperCRUNCH consists of a set of modules written in Python that function as stand-alone command-line scripts. These modules are available in the [supercrunch-scripts](https://github.com/dportik/SuperCRUNCH/tree/master/supercrunch-scripts) folder. They can be downloaded and executed independently without the need to install SuperCRUNCH as a Python package or library, making them easy to use and edit. The scripts function independently, and do not require being contained or used in the same directory. There are several external dependencies that should be installed prior to use of SuperCRUNCH if you plan to use all the available modules, including:

+ [**Biopython**](https://biopython.org/)
+ [**NCBI-BLAST+**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ [**CD-HIT-EST**](http://weizhongli-lab.org/cd-hit/)
+ [**MAFFT**](https://mafft.cbrc.jp/alignment/software/)
+ [**Muscle**](https://www.drive5.com/muscle/)
+ [**Clustal-O**](http://www.clustal.org/omega/)
+ [**MACSE**](https://bioweb.supagro.inra.fr/macse/)
+ [**trimAl**](http://trimal.cgenomics.org/)

Helpful installation instructions for these dependencies can be found on the [**Installation Instructions**](https://github.com/dportik/SuperCRUNCH/wiki/Installation-Instructions) wiki. Please note that some modules do not require any dependencies. If you plan to use only a subset of modules you can quickly check which modules require dependencies [here](https://github.com/dportik/SuperCRUNCH/wiki/Installation-Instructions#module-dependencies-list). 

SuperCRUNCH scripts can be run using Mac OSX (10.10+) and Linux, and can also work with Windows using a program like Cygwin. 


## Complete Instructions for Analyses

SuperCRUNCH has extensive documentation which can be accessed through the wiki tab at the top of the page. An overview of the components of SuperCRUNCH can be found on the [**Analysis Overview**](https://github.com/dportik/SuperCRUNCH/wiki/1:-Analysis-Overview) page. This section outlines all major topics and navigates to detailed instructions for each step, including usage for all modules, proposed workflows, and common issues. Quick links to wiki pages for the major components of SuperCRUNCH are also provided below:

+ [**Starting Sequences**](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Sequences)
+ [**Assessing Taxonomy**](https://github.com/dportik/SuperCRUNCH/wiki/3:-Assessing-Taxonomy)
+ [**Parsing Loci**](https://github.com/dportik/SuperCRUNCH/wiki/4:-Parsing-Loci)
+ [**Similarity Searches and Filtering**](https://github.com/dportik/SuperCRUNCH/wiki/5:-Similarity-Searches-and-Filtering)
+ [**Sequence Selection**](https://github.com/dportik/SuperCRUNCH/wiki/6:-Sequence-Selection)
+ [**Multiple Sequence Alignment**](https://github.com/dportik/SuperCRUNCH/wiki/7:-Multiple-Sequence-Alignment)
+ [**Post-Alignment Tasks**](https://github.com/dportik/SuperCRUNCH/wiki/8:-Post-Alignment-Tasks)


## Tutorials and Examples

New tutorials that use SuperCRUNCH v1.2 are under construction, please check back soon! 

## Citation

SuperCRUNCH is currently described in a pre-print available on BioRxiv:

+ Portik, D.M., and J.J. Wiens. (2019) SuperCRUNCH: A toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. BioRxiv, https://doi.org/10.1101/538728.

SuperCRUNCH is also in peer-review, and we hope to have a formal publication soon. 

If you use SuperCRUNCH for your research, please cite the above BioRxiv publication for now.

## License

GNU Lesser General Public License v3.0

## Contact

SuperCRUNCH is written and maintained by Daniel Portik (daniel.portik@gmail.com)
