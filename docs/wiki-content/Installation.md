## Installing Dependencies

**SuperCRUNCH** requires several external dependencies for full functionality, including:

+ [Biopython](https://biopython.org/)
+ [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ [CD-HIT-EST](http://weizhongli-lab.org/cd-hit/)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/)
+ [Muscle](https://www.drive5.com/muscle/)
+ [Clustal-O](http://www.clustal.org/omega/)
+ [MACSE](https://bioweb.supagro.inra.fr/macse/)
+ [trimAl](http://trimal.cgenomics.org/)

### ***Python Packages***

**SuperCRUNCH** is written in Python (2.7). It requires having the [Numpy](http://www.numpy.org/) and [Biopython]((https://biopython.org/)) python packages installed. Numpy is a package for scientific computing, and it is generally included with the Python library in newer platforms. Numpy can also be installed as part of the SciPy package collection and instructions for doing so can be found [here](http://www.numpy.org/). Biopython is used to read and write sequence data files and perform various functions. Biopython can be installed using pip or installed from source, and instructions for doing so can be found [here](https://biopython.org/wiki/Download). **SuperCRUNCH** has been tested for BioPython v1.7.2 and 1.7.3.

### ***External Dependencies***

**SuperCRUNCH** also relies on several additional programs, including [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (specifically the *blastn* and *makeblastdb* programs), [CD-HIT-EST](http://weizhongli-lab.org/cd-hit/), [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Muscle](https://www.drive5.com/muscle/), [Clustal-O](http://www.clustal.org/omega/), [MACSE](https://bioweb.supagro.inra.fr/macse/) (v2), and [trimAl](http://trimal.cgenomics.org/). Download and installation instructions can be found using the above links for each program. With the exception of MACSE, which is a jar file (java executable), all other programs need to be in PATH to work properly for **SuperCRUNCH**. This can be accomplished by putting all the executables in a directory and setting it in PATH (quick guide [here](http://osxdaily.com/2014/08/14/add-new-path-to-path-command-line/)). 

**NOTE:** I encountered serious problems when I installed and used the openmp version of cd-hit-est for automated clustering in *Cluster_Blast_Extract.py* - the clusters created were inconsistent and I could not replicate results across identical runs. This issue was resolved using the regular version of cd-hit-est, which was compiled using 'make openmp=no'. I *strongly* recommend doing the same, and following the instructions on the website to accomplish this.

After installation, please check to make sure that the relevant executables are in PATH (type their name on the command line and they should run) and that they have the following names which are case-sensitive (if not, relabel them accordingly):

+ **NCBI-BLAST+**: *blastn*, *makeblastdb*, tested for version 2.7.1+
+ **CD-HIT-EST**: *cd-hit-est*, tested for version 4.7 (non-openmp)
+ **MAFFT**: *mafft*, tested for version 7.407
+ **Muscle**: *muscle*, tested for version 3.8.31
+ **Clustal-O**: *clustalo*, tested for version 1.2.3
+ **trimAl**: *trimal*, tested for version 1.4.22

This should allow everything to run properly. For quick reference, here is a list of the currently available **SuperCRUNCH** modules and their dependencies:

+ **Adjust_Direction.py**: *Biopython*, *mafft*
+ **Align.py**: *Biopython*, *mafft*, *muscle*, *clustalo*, *macse*
+ **Cluster_Blast_Extract.py**: *Biopython*, *cd-hit-est*, *blastn*, *makeblastdb*
+ **Coding_Translation_Tests.py**: *Biopython*
+ **Concatenation.py**: None!
+ **Contamination_Filter.py**: *Biopython*, *blastn*, *makeblastdb*
+ **Fasta_Convert.py**: *Biopython*
+ **Filter_Seqs_and_Species.py**: *Biopython*
+ **Make_Acc_Table.py**: None!
+ **Parse_Loci.py**: *Biopython*, *Numpy*
+ **Reference_Blast_Extract.py**: *Biopython*, *blastn*, *makeblastdb*
+ **Relabel_Fasta.py**: None!
+ **Remove_Duplicate_Accessions.py**: *Biopython*
+ **Rename_Merge.py**: *Biopython*, *Numpy*
+ **Taxa_Assessment.py**: *Biopython*, *Numpy*
+ **Trim_Alignments.py**: *Biopython*, *trimal*

-----------

*Last updated: January 2019*