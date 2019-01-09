![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Overview

**SuperCRUNCH** is a python program used to create, filter and manipulate large phylogenetic data sets based on nucleotide data downloaded directly from NCBI. It can be used to search for specified sets of taxa and loci to create population or species-level molecular data sets, including phylogenetic supermatrices. It offers a variety of transparent and flexible options for orthology detection, sequence selection, alignment, and file manipulation. The main **SuperCRUNCH** workflow allows rapid construction of phylogenetic supermatrices from downloaded sequence data, greatly simplifying the process of updating and improving large phylogenies with newly available molecular data. Because **SuperCRUNCH** streamlines tasks like alignment, trimming, and concatenation, it can also be used to efficiently process sequence capture data sets. 

## Installation

**SuperCRUNCH** consists of several Python modules that can be used as command line scripts, and therefore does not require any installation. The scripts function independently, and do not require being used in the same directory. There are several external dependencies that should be installed prior to use of **SuperCRUNCH**, including:

+ [Biopython](https://biopython.org/)
+ [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ [CD-HIT-EST](http://weizhongli-lab.org/cd-hit/)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/)
+ [Muscle](https://www.drive5.com/muscle/)
+ [Clustal-O](http://www.clustal.org/omega/)
+ [MACSE](https://bioweb.supagro.inra.fr/macse/)
+ [trimAl](http://trimal.cgenomics.org/)

Helpful installation instructions for these dependencies can be found on the wiki page [here](to fill in).


## Tutorials and Example Analyses

Complete instructions for running **SuperCRUNCH** analyses are available on the wiki page [here](to fill in), including detailed usage for all modules, proposed workflows, and common issues.

Several complete example analyses are available two unrelated clades:

+ ***Iguania*** - an infraorder of squamate reptiles
+ ***Dipsacales*** - an order of flowering plants 

The complete set of material for these analyses, including all input files, command instructions, and output files, is freely available [here](https://osf.io/bpt94/) on the **SuperCRUNCH** Open Science Framework project page.

## Citation

**SuperCRUNCH** was introduced in the following publication:


If you use **SuperCRUNCH** for your research, please cite this publication.

## License

GNU Lesser General Public License v3.0

## Contact

Daniel Portik, PhD

daniel.portik@gmail.com
