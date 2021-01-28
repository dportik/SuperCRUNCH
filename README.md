![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Overview

**SuperCRUNCH** is a python toolkit for creating and working with phylogenetic datasets. SuperCRUNCH can be run using any set of sequence data, as long as sequences are in fasta format with standard naming conventions (described [here](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Sequences)). 

**SuperCRUNCH** can be used to:
+ Construct de novo supermatrices from GenBank and/or local sequence data.
+ Construct phylogeographic datasets from GenBank and/or local sequence data, with the ability to detect voucher codes.
+ Parse a large fasta file of GenBank and/or local sequences into gene-specific fasta files, based on a list of genes and list of species names.
+ Perform similarity filtering for all sequences of a gene to remove mis-identified or highly divergent sequences.
+ Extract specific genes from a set of mitochondrial genomes.
+ Use a set of reference sequences to trim all input sequences to match the reference region. 
+ Adjust reading frames for coding sequences, and identify problematic coding sequences.
+ Adjust sequence directions to prepare for alignment.
+ Automate multiple sequence alignment using Clustal-O, MAFFT, Muscle, and MACSE.
+ Generate accession tables for supermatrices, for all genes and species included.
+ Relabel GenBank sequences using species names, accession numbers, and/or voucher codes.
+ Trim alignments using trimAl or a custom trimming routine
+ Convert fasta format to phylip and nexus
+ Concatenate any number of fasta or phylip alignment files to create concatenated alignments
+ And many other tasks!


A visual overview of the major steps in SuperCRUNCH is shown below:

![SuperCrunch workflow](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Figure-1.jpg)

SuperCRUNCH is highly modular and analyses do not require running the full pipeline. There are many useful tools available for manipulating your phylogenetic and phylogeographic datasets.

## Citation 

SuperCRUNCH is described in more detail in the following publication:

+ Portik, D.M., and J.J. Wiens. (2020) SuperCRUNCH: A bioinformatics toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. Methods in Ecology and Evolution, 11: 763-772. https://doi.org/10.1111/2041-210X.13392

The published article is available [**here**](https://github.com/dportik/SuperCRUNCH/tree/master/docs/publication). A pre-print was made available on BioRxiv prior to publication ([**here**](https://www.biorxiv.org/content/10.1101/538728v3)).


## Installation

There are several dependencies required to run SuperCRUNCH, including Python packages ([**BioPython**](https://biopython.org/) and [numpy](https://numpy.org/)), as well as external tools ([**NCBI-BLAST+**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [**CD-HIT-EST**](http://weizhongli-lab.org/cd-hit/), [**MAFFT**](https://mafft.cbrc.jp/alignment/software/), [**Muscle**](https://www.drive5.com/muscle/), [**Clustal-O**](http://www.clustal.org/omega/), [**MACSE**](https://bioweb.supagro.inra.fr/macse/), and [**trimAl**](http://trimal.cgenomics.org/)). 

Installation of these requirements is fast and easy using `conda`. The [supercrunch-conda-env.yml](https://github.com/dportik/SuperCRUNCH/blob/master/supercrunch-conda-env.yml) file can be used to create the correct conda environment:

```
conda create -f supercrunch-conda-env.yml
```

The resulting conda environment can then be activated using:

```
conda activate supercrunch
```

You can then run all SuperCRUNCH modules in this environment. 

**NOTE** - this will install all requirements except for `MACSE`, which is a jar file that must be downloaded from [**here**](https://bioweb.supagro.inra.fr/macse/index.php?menu=releases) (get V2.05). This software is only required for the MACSE option in the alignment module (`Align.py`).

For non-conda installation of these packages, please see the [**Installation Instructions**](https://github.com/dportik/SuperCRUNCH/wiki/Installation-Instructions) wiki. 

SuperCRUNCH itself consists of a set of modules written in Python (compatible with 2.7 and 3.7) that function as stand-alone command-line scripts. These modules are available in the [supercrunch-scripts](https://github.com/dportik/SuperCRUNCH/tree/master/supercrunch-scripts) folder. They can be downloaded and executed independently without the need to install SuperCRUNCH as a Python package or library, making them easy to use and edit. The scripts function independently, and do not require being contained or used in the same directory. SuperCRUNCH scripts can be run using Mac OSX (10.10+) and Linux, and can also work with Windows using a program like Cygwin. 

## Version

The current release of **SuperCRUNCH** is [**v1.3.0**](https://github.com/dportik/SuperCRUNCH/releases). Please see below for important changes.
 
#### Changes in v1.3.0:
  - Added a `conda` environment recipe for SuperCRUNCH, allowing easy installation of all requirements except MACSE.
  - `Parse_Loci.py`: Added new feature that allows a term to be added to the loci search terms that will exclude a record if a match is found. For example, adding the negative term `pseudogene` will exclude all records containing that word, even if they match the other abbreviation or description terms. This requires a four-column search terms file, where the fourth column is the negative term (`N/A` in this column indicates no negative term should be used). This module was made backwards-compatible with the three-column search terms file - if a fourth column is not present the `N/A` is automatically generated.
  - `Filter_Seqs_and_Species.py`: Added `--accessions_include` flag. This points to a text file of accession numbers (one per line). When used with the `--seq_selection oneseq` option, if an accession included in the list is found in the available seqs for a taxon and gene, it must be selected. This is not just an "allowed list", this list will override other settings for selection such as length. Also added the `--accessions_exclude` flag, which points to a text file of accession numbers (one per line). These accessions will NEVER be selected - they are removed from all searches. This is the equivalent of including a "blocked list".
  - `Taxa_Assessment.py`: Altered SQL search query for "unmatched" taxa to avoid sql variable limit maximum issue. Also, now invokes the `SeqIO.index_db()` method for sequence files >5GB, rather than using `SeqIO.index()` method, which is much more memory efficient for big data. The `SeqIO.index_db()` method is already used in `Parse_Loci.py`.
  - `Cluster_Blast_Extract.py`: Added feature to remove problematic long sequences if they somehow end up in the main cluster of sequences for a gene. The new filter removes all seqs that are 1.3x the length of the 95th percentile of all lengths.
  - Added a new `Remove_Long_Accessions.py` module, which can filter a downloaded GenBank fasta file to remove extremely long sequences (>150kb). This will eliminate whole genome sequencing records, which are not useful for SuperCRUNCH.
  - Updated recognition for file extensions produced by updated blastn tools (`.ndb`, `.not`, `.ntf`, `.nto`).

For complete version history please see the [releases](https://github.com/dportik/SuperCRUNCH/releases).


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

Several tutorials were made available as part of the original SuperCRUNCH publication, which cover the full range of analyses available. These tutorials can be found on the [OSF SuperCRUNCH project page](https://osf.io/bpt94/). Each tutorial includes all data and instructions necessary to replicate the analysis. An overview of the tutorials available can be found on the [Tutorials wiki](https://github.com/dportik/SuperCRUNCH/wiki/Tutorials).

## Reporting Issues, Getting Help, and Providing Suggestions

For main analysis issues and/or bugs, please create an issue on github [here](https://github.com/dportik/SuperCRUNCH/issues). Make sure you include the details of your analysis (inputs, outputs, commands) to assist the troubleshooting. The [**SuperCRUNCH user group**](http://groups.google.com/group/supercrunch-users) can also be used to create a post.


## License

GNU Lesser General Public License v3.0

## Contact

SuperCRUNCH is written and maintained by Daniel Portik (daniel.portik@gmail.com)
