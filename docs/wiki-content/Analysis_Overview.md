![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

# Analysis Overview

**SuperCRUNCH** is a python toolkit for extracting, filtering, and manipulating nucleotide data. It is modular in design, and there are many ways to perform analyses to meet specific research goals. This page provides an overview of the major topics involved in a full **SuperCRUNCH** analysis, recommendations for performing analyses, and example workflows. 

## Summary of SuperCRUNCH Components

Below, a summary of the main topics and scripts is provided in the approximate order in which they would be used. Complete instructions for performing each step are provided in the individual wiki pages for each topic. These can be reached by clicking on the relevant link below or to the right on the wiki toolbar. **Helpful information can always be displayed on the command line by running any script using the -h flag.** 

### [Starting Materials](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Materials)

+ Obtaining Sequence Data
    + Using Custom Sequences
    + Remove_Duplicate_Accessions.py
+ Obtaining Taxon Names Lists
    + Getting Taxa From Fasta Files
+ Obtaining Loci Search Terms
    + Searching for UCE loci

### [Taxon Filtering and Locus Parsing](https://github.com/dportik/SuperCRUNCH/wiki/3:-Taxon-Filtering-and-Locus-Parsing)

+ Taxa_Assessment.py
+ Rename_Merge.py
+ Parse_Loci.py

### [Orthology Filtering](https://github.com/dportik/SuperCRUNCH/wiki/4:-Orthology-Filtering)

+ Cluster_Blast_Extract.py
+ Reference_Blast_Extract.py
+ Contamination_Filter.py

### [Sequence Quality Filtering and Selection](https://github.com/dportik/SuperCRUNCH/wiki/5:-Sequence-Quality-Filtering-and-Selection)

+ Filter_Seqs_and_Species.py
+ Make_Acc_Table.py
+ Infer_Supermatrix_Combinations.py

### [Sequence Alignment](https://github.com/dportik/SuperCRUNCH/wiki/6:-Sequence-Alignment)

+ Adjust_Direction.py
+ Coding_Translation_Tests.py
+ Align.py

### [Post-Alignment Tasks](https://github.com/dportik/SuperCRUNCH/wiki/7:-Post-Alignment-Tasks)

+ Relabel_Fasta.py
+ Trim_Alignments.py
+ Fasta_Convert.py
+ Concatenation.py

A typical **SuperCRUNCH** run includes executing a majority of these steps. However, this workflow is extremely flexible and can be tailored to achieve a variety of goals.

## Analysis Recommendations

There are many ways to perform analyses using SuperCRUNCH, and you may find yourself using all the modules or perhaps just one. Regardless of your purpose for using SuperCRUNCH, there are a few helpful suggestions to keep in mind when running your analyses:

+ **Organization.** Consider how you will organize your starting input files, and how you will manage files generated across steps. 
+ **File Management.** It is good practice to run a module, identify the key output files, and then move or copy/paste these files into a new directory before running any following steps. Numbering the new directories is a good way to keep track of the progression of your analysis, similar to what is illustrated in the OSF [example analyses](https://osf.io/bpt94/). For some modules you are required to specify the output directory, but for other modules the output directories are written automatically. Make sure you look at the summary of outputs  provided in the detailed wiki instructions for each module. 
+ **Input File Formats and Example Files.** If you are not planning on using the full pipeline, please pay careful attention to how you format your input files and also how you name them. Refer to the detailed wiki instructions for each module. For examples of various files, see the [example data folder](https://github.com/dportik/SuperCRUNCH/tree/master/data). The OSF [example analyses](https://osf.io/bpt94/) will have all the input and output files from every step of the analysis, and can be referred to for examples of all types of input files.
+ **Run Errors.** If a module crashes or ends with an error, you should delete any intermediate/output files created before attempting to run the module again. By default, SuperCRUNCH will append to files rather than overwrite them, which can introduce unwanted errors in the output if a script is run repeatedly.
+ **Error Reporting.** If you encounter any bugs or major problems while running SuperCRUNCH, please post them to the [issues](https://github.com/dportik/SuperCRUNCH/issues) page on github so they can be addressed. 

## Example Workflows

Examples of various workflows are provided here. The complete workflows assume that all starting materials have been gathered previously (starting sequence data, taxon names, loci list). An example of a more creative use of SuperCRUNCH is provided at the bottom.

### ***de novo Supermatrix***

This analysis workflow will result in a single representative sequence per taxon per locus, which is the traditional supermatrix used for inferring species-level phylogenies. It is strongly recommended that the output files of each step be moved or copy/pasted into a new directory before performing the subsequent step.

+ Perform an initial search for taxon names in the starting sequence set using `Taxa_Assessment.py`. 
+ _Optional:_ Correct taxon names using `Rename_Merge.py`.
+ Create individual fasta files for all loci using `Parse_Loci.py`.
+ Perform orthology-filtering for sequences within the locus-specific fasta files using `Cluster_Blast_Extract.py` and/or `Reference_Blast_Extract.py`, based on the type of records contained within the fasta files. 
+ _Optional:_ Screen sequences for possible contamination using `Contamination_Filter.py`. For example, in a set of reptile mtDNA sequences find and remove any human contamination.
+ Select the best representative sequence per taxon per locus using `Filter_Seqs_and_Species.py`, based on desired qualities (length, translatability, random).
+ Create a table of the accession numbers from all the selected sequences across loci using `Make_Acc_Table.py`.
+ _Optional (but very fun):_ Infer the total number of possible supermatrix combinations given the total filtered sequence set using `Infer_Supermatrix_Combinations.py`.
+ Adjust incorrect sequence directions within the filtered fasta files prior to alignment using `Adjust_Direction.py`.
+ Align sequences using one or all of the available aligners in `Align.py`. For translation alignments, running `Coding_Translation_Tests.py` is necessary prior to aligning.
+ Relabel sequences from the original NCBI description lines to taxon names using `Relabel_Fasta.py` with the `-r species` flag.
+ _Optional:_ Trim alignments using `Trim_Alignments.py `.
+ Concatenate the sequences of all loci to create a supermatrix using `Concatenation.py`.
+ Use final supermatrix for phylogenetic analyses.


### ***Population-level data***

This analysis workflow will result in multiple sequences per taxon per locus, which can be used for phylogeographic or population genetic analyses. It is strongly recommended that the output files of each step be moved or copy/pasted into a new directory before performing the subsequent step.

+ Perform an initial search for taxon names in the starting sequence set using `Taxa_Assessment.py`. 
+ _Optional:_ Correct taxon names using `Rename_Merge.py`.
+ Create individual fasta files for all loci using `Parse_Loci.py`.
+ Perform orthology-filtering for sequences within the locus-specific fasta files using `Cluster_Blast_Extract.py` and/or `Reference_Blast_Extract.py`, based on the type of records contained within the fasta files. 
+ _Optional:_ Screen sequences for possible contamination using `Contamination_Filter.py`. For example, in a set of reptile mtDNA sequences find and remove any human contamination.
+ Select all sequences per taxon per locus that pass the minimum length filter using `Filter_Seqs_and_Species.py` with the `--allseqs` flag.
+ Adjust incorrect sequence directions within the filtered fasta files prior to alignment using `Adjust_Direction.py`.
+ Align sequences using one or all of the available aligners in `Align.py`. For translation alignments, running `Coding_Translation_Tests.py` is necessary prior to aligning.
+ Relabel sequences from the original NCBI description lines to taxon names plus accessions using `Relabel_Fasta.py` with the `-r species_acc` flag.
+ _Optional:_ Trim alignments using `Trim_Alignments.py `.
+ Use loci for gene tree construction or various population genetics analyses.

### ***Preparing protein-coding loci for NCBI submission to BankIt***

To submit protein-coding sequences to NCBI, they should ideally all begin in the same reading frame. Any internal stop codons are likely to cause problems with submissions, so they should be identified and corrected beforehand. The following is an example of how to use SuperCRUNCH to process unaligned fasta files in preparation for submission.

+ Use `Coding_Translation_Tests.py` to identify a correct reading frame of each sequence within the unaligned fasta file.
+ The output file labeled `[Input Fasta name]_Passed.fasta` will contain all sequences that passed translation (no internal stop codons), which have each been adjusted to reading frame 1.
+ The output file labeled `[Input Fasta name]_Failed.fasta` will contain all sequences that failed translation (contain at least one internal stop codon), which will need to be inspected for errors.
+ If you believe any sequences in the `[Input Fasta name]_Failed.fasta` are a result of sequencing errors, you can attempt to auto-correct these by running `Align.py` using the paired alignment feature of MACSE. The resulting alignment will contain sequences with error-free reading frames with all sequences beginning in frame 1. The alignment can be stripped and tested again with `Coding_Translation_Tests.py` to ensure all sequences begin in frame 1 and do not contain any internal stop codons.
+ The `Coding_Translation_Tests.py` sometimes introduces N's at the end of sequences to complete the final codon. These should be trimmed prior to submission. 

