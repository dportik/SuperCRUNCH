![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Overview

**SuperCRUNCH** is a python toolkit for creating and working with phylogenetic datasets. SuperCRUNCH can be run using any set of sequence data, as long as sequences are in fasta format with standard naming conventions (described [here](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Materials)). 

SuperCRUNCH can be used to process sequences downloaded directly from GenBank/NCBI, local sequence data (e.g. sequences not downloaded from GenBank, such as unpublished data), or a combination of both. The sequence data are first parsed into gene-specific fasta files using targeted searches guided by lists of taxa and locus names. For each resulting gene, sequences can be filtered with similarity searches using automated methods or based on user-supplied reference sequences. SuperCRUNCH offers the option to select a best representative sequence for each taxon, or to retain all filtered sequences for each taxon. These options allow the user to generate interspecific supermatrix datasets (one sequence per taxon per locus) or population-level datasets (multiple sequences per taxon per locus). SuperCRUNCH offers important pre-alignment steps (adjust sequence directions, adjust reading frames) and several options for multiple sequence alignment (Clustal-O, MAFFT, Muscle, MACSE) and alignment trimming (using trimAl). Finally, SuperCRUNCH can be used for rapid file format conversion and concatenation. 

SuperCRUNCH is scalable and can be used to assemble a variety of datasets, ranging from small population-level datasets (one taxon, one gene) to large phylogenomic datasets with thousands of loci (such as UCEs or other sequence capture datasets). SuperCRUNCH was intended to be transparent, objective and repeatable, and provides meaningful output at every step to help guide user decisions. In addition, it is modular in design and various components of SuperCRUNCH can be easily incorporated into custom bioinformatics workflows.

A full overview of **SuperCRUNCH** is described in the following pre-print article:

+ Portik, D.M., and J.J. Wiens. (2019) SuperCRUNCH: A toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. BioRxiv, https://doi.org/10.1101/538728.


## Version

The current release of **SuperCRUNCH** is [v1.1](https://github.com/dportik/SuperCRUNCH/releases). For version history please see the [change log file](https://github.com/dportik/SuperCRUNCH/tree/master/CHANGELOG.md).


## Installation

**SuperCRUNCH** consists of a set of modules written in Python (2.7) that function as stand-alone command-line scripts. These modules are available in the [supercrunch-scripts](https://github.com/dportik/SuperCRUNCH/tree/master/supercrunch-scripts) folder. They can be downloaded and executed independently without the need to install SuperCRUNCH as a Python package or library, making them easy to use and edit. The scripts function independently, and do not require being contained or used in the same directory. There are several external dependencies that should be installed prior to use of SuperCRUNCH if you plan to use all the available modules, including:

+ [Biopython](https://biopython.org/)
+ [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ [CD-HIT-EST](http://weizhongli-lab.org/cd-hit/)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/)
+ [Muscle](https://www.drive5.com/muscle/)
+ [Clustal-O](http://www.clustal.org/omega/)
+ [MACSE](https://bioweb.supagro.inra.fr/macse/)
+ [trimAl](http://trimal.cgenomics.org/)

Helpful installation instructions for these dependencies can be found on the [wiki installation page](https://github.com/dportik/SuperCRUNCH/wiki/Installation-Instructions).
Please note that some modules do not require any dependencies. If you plan to use only a subset of modules you can quickly check which modules require dependencies [here](https://github.com/dportik/SuperCRUNCH/wiki/Installation-Instructions). 

SuperCRUNCH scripts can be run using Mac OSX (10.10+) and Linux, and can also work with Windows if using Cygwin. 


## Complete Instructions for Analyses

An overview of the components of **SuperCRUNCH** can be found on the [wiki overview page](https://github.com/dportik/SuperCRUNCH/wiki/1:-Analysis-Overview). This page outlines all major topics and navigates to detailed instructions for each step, including usage for all modules, proposed workflows, and common issues.

## Tutorials and Examples

Several example analyses associated with the pre-print are available on the [**SuperCRUNCH** project page](https://osf.io/bpt94/) on the Open Science Framework, including all input files, output files, and complete wiki tutorials. These analyses include:

+ **Iguania Supermatrix:** SuperCRUNCH analysis of Iguania using sequence data downloaded directly from NCBI. [Direct link](https://osf.io/vsu5k/).
+ **UCE Supermatrix:** SuperCRUNCH analysis of UCE loci available for the frog genus *Kaloula*, using sequence data downloaded directly from NCBI. [Direct link](https://osf.io/5rey2/).
+ **Population Dataset (*Uma*):** SuperCRUNCH analysis to retrieve population-level data for the lizard genus *Uma*, using sequence data downloaded directly from NCBI. [Direct link](https://osf.io/49wz3/).
+ **Population Dataset (*Callisaurus*):** SuperCRUNCH analysis to retrieve population-level data for the lizard genus *Callisaurus*, using sequence data downloaded directly from NCBI. [Direct link](https://osf.io/6bwhf/).

On the OSF project page, there are also several analyses that compare the performance of **SuperCRUNCH** to other programs. 

The analyses available on the OSF page will have examples of all the input files required to run various steps. Additionally, several example input files are also provided in the [data folder](https://github.com/dportik/SuperCRUNCH/tree/master/data). These include several locus search terms files (including the 5k UCE set) and multiple reference sequence sets. 

**SuperCRUNCH** was presented at the [**Trees in the Desert workshop**](http://db.herbarium.arizona.edu/phlora/TID/WorkshopPage.html), which was held to discuss challenges related to ultra-large phylogenetic trees. Several materials from the workshop, including a [presentation](https://github.com/dportik/SuperCRUNCH/tree/master/docs/Trees-in-the-desert-workshop/Portik_SuperCRUNCH_presentation.pdf) and short tutorial (with input data) are available [here](https://github.com/dportik/SuperCRUNCH/tree/master/docs/Trees-in-the-desert-workshop/). 

## Citation

**SuperCRUNCH** is currently described in a pre-print available on BioRxiv:

+ Portik, D.M., and J.J. Wiens. (2019) SuperCRUNCH: A toolkit for creating and manipulating supermatrices and other large phylogenetic datasets. BioRxiv, https://doi.org/10.1101/538728.

SuperCRUNCH is also in peer-review, and we hope to have a formal publication soon. 

If you use SuperCRUNCH for your research, please cite the above BioRxiv publication (for now).

## License

GNU Lesser General Public License v3.0

## Contact

SuperCRUNCH is written and maintained by Daniel Portik (daniel.portik@gmail.com)
