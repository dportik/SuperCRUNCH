![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Overview

**SuperCRUNCH** is a python toolkit for creating and working with phylogenetic datasets. SuperCRUNCH can be run using any set of sequence data, as long as sequences are in fasta format with standard naming conventions (described [here](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Sequences)). 

SuperCRUNCH can be used to process sequences downloaded directly from GenBank/NCBI, local sequence data (e.g. sequences not downloaded from GenBank, such as unpublished data), or a combination of both. The sequence data are first parsed into gene-specific fasta files using targeted searches guided by lists of taxon and locus names. For each resulting gene, sequences can be filtered with similarity searches using automated methods or based on user-supplied reference sequences. SuperCRUNCH offers the option to select a best representative sequence for each taxon, or to retain all filtered sequences for each taxon. These options allow the user to generate species-level supermatrix datasets (one sequence per species per locus) or population-level datasets (multiple sequences per species per locus). In addition, SuperCRUNCH can identify voucher codes present in sequence records and link samples in phylogeographic datasets through correct labeling. SuperCRUNCH offers important pre-alignment steps (adjust sequence directions, adjust reading frames), several options for sequence alignment (Clustal-O, MAFFT, Muscle, MACSE), and multiple options for alignment trimming. Finally, SuperCRUNCH can be used for rapid file format conversion and concatenation. A visual overview of the major steps in SuperCRUNCH is shown below:

![SuperCrunch workflow](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Figure-1.jpg)

SuperCRUNCH is scalable and can be used to assemble a variety of datasets, ranging from small population-level datasets (one taxon, one gene) to large phylogenomic datasets with thousands of loci (such as UCEs or other sequence capture datasets). SuperCRUNCH was intended to be transparent, objective and repeatable, and provides meaningful output at every step to help guide user decisions. In addition, it is modular in design and various components of SuperCRUNCH can be easily incorporated into other custom bioinformatics workflows.

SuperCRUNCH is described in more detail in the following publication:

+ Portik, D.M., and J.J. Wiens. (2020) SuperCRUNCH: A bioinformatics toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. Methods in Ecology and Evolution, 11: 763-772. https://doi.org/10.1111/2041-210X.13392

The article is available [**here**](https://github.com/dportik/SuperCRUNCH/tree/master/docs/publication), and the pre-print is available on BioRxiv [**here**](https://www.biorxiv.org/content/10.1101/538728v3).


## Installation

SuperCRUNCH consists of a set of modules written in Python (compatible with 2.7 and 3.7) that function as stand-alone command-line scripts. These modules are available in the [supercrunch-scripts](https://github.com/dportik/SuperCRUNCH/tree/master/supercrunch-scripts) folder. They can be downloaded and executed independently without the need to install SuperCRUNCH as a Python package or library, making them easy to use and edit. The scripts function independently, and do not require being contained or used in the same directory. There are several external dependencies that should be installed prior to use of SuperCRUNCH if you plan to use all the available modules, including:

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


## Got a question, need some help, or have a suggestion?

Please head to the [**SuperCRUNCH user group**](http://groups.google.com/group/supercrunch-users) and create a post. The user group can be found here: http://groups.google.com/group/supercrunch-users. 



## Tutorials and Examples

There are currently nine tutorials, which cover the full range of analyses run in the SuperCRUNCH publication. These tutorials can be found on the [OSF SuperCRUNCH project page](https://osf.io/bpt94/). Each tutorial includes all data and instructions necessary to replicate the analysis. The tutorials include:

- [**Iguania-Fast**](https://osf.io/x5hrm/): **Construct a very large species-level supermatrix for Iguania using "fast" settings.**
    - Produces supermatrix with 60 loci, 1,399 species, 12,978 sequences.
    - Takes ~1.5 hours to run and requires 15GB of space. 
    - Difficulty: Medium.
- [**Iguania-Thorough**](https://osf.io/9gs32/): **Construct a very large species-level supermatrix for Iguania using "thorough" settings.** 
    - Produces supermatrix with 61 loci, 1,426 species, 13,307 sequences.
    - Takes ~13 hours to run and requires 30GB of space. 
    - Difficulty: Hard.
- [**Kaloula-Vouchered**](https://osf.io/zxnq8/): **Construct a vouchered population-level UCE phylogenomic supermatrix.**
    - Produces supermatrix with 1,784 loci, 18 vouchered samples, 28,790 sequences.
    - Takes ~20 minutes to run and requires 350MB of space. 
    - Difficulty: Easy.
- [**Kaloula-Species**](https://osf.io/crzp5/): **Construct a species-level UCE phylogenomic supermatrix, which is a variation of the above tutorial.** 
    - Produces supermatrix with 1,784 loci, 14 species, 22,717 sequences.
    - Takes ~20 minutes to run and requires 350MB of space. 
    - Difficulty: Easy.
- [**Trachylepis-Vouchered**](https://osf.io/bgc5z/): **Reconstruct a vouchered multi-locus phylogeography dataset from a published study.** 
    - Produces dataset with 4 loci, 108 vouchered samples, 400 sequences.
    - Takes ~1 minute to run and requires 10MB of space. 
    - Difficulty: Easy.
- [**Trachylepis-Species**](https://osf.io/umswn/): **Construct a species-level supermatrix from a published phylogeography dataset.** 
    - Produces supermatrix with 4 loci, 7 species, 26 sequences.
    - Takes ~1 minute to run and requires 10MB of space. 
    - Difficulty: Easy.
- [**Callisaurus-Population**](https://osf.io/7gujb/): **Construct a de novo multi-locus population-level dataset from multiple published sources.** 
    - Produces dataset with 7 loci and 265 sequences. Samples are NOT linked by voucher codes.
    - Takes ~1 minute to run and requires 150MB of space. 
    - Difficulty: Easy.
- [**Uma-Population**](https://osf.io/e28tu/): **Construct a de novo multi-locus population-level dataset from multiple published sources.** 
    - Produces dataset with 5 loci and 234 sequences. Samples are NOT linked by voucher codes.
    - Takes ~1 minute to run and requires 75MB of space. 
    - Difficulty: Easy.
- [**Hyperoliid-Outgroup**](https://osf.io/q9nyx/): **Incorporate published outgroup sequences from GenBank with unpublished ingroup sequences to create a custom supermatrix.** 
    - Produces supermatrix with 6 loci, 365 samples, 1,724 sequences.
    - Takes ~5 minutes to run and requires 100MB of space. 
    - Difficulty: Medium.


## Citation

SuperCRUNCH was published as:

+ Portik, D.M., and J.J. Wiens. (2020) SuperCRUNCH: A bioinformatics toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. Methods in Ecology and Evolution, 11: 763-772. https://doi.org/10.1111/2041-210X.13392


## Version

The current release of **SuperCRUNCH** is [**v1.2.1**](https://github.com/dportik/SuperCRUNCH/releases). 

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


## License

GNU Lesser General Public License v3.0

## Contact

SuperCRUNCH is written and maintained by Daniel Portik (daniel.portik@gmail.com)
