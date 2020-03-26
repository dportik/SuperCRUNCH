## Tutorial for Analysis: Iguania-Thorough

This wiki is intended to allow you to replicate the Iguania-Thorough analysis reported in the Supplemental Materials. All of the input and output files of the completed analysis are provided in the OSF data folders.  

The wiki is presented as a tutorial. However, you can always scroll through the sections of this wiki to see how each module was called, and which files were used and produced at each step. One of the main goals is to provide you with all the input files required for running the different modules of SuperCRUNCH, and to demonstrate the use of various options along the way. This page is therefore intended to function as a tutorial and a reference resource.

**Please note that this is the most complicated analysis that we ran. Running all steps in the tutorial will take at least 13 hours, and requires ~50 GB of free space. If you'd like a less complicated or smaller sized tutorial, you should look at any other analysis besides this one!**

---------------


## **Preparing for the Tutorial**

Before you can run this tutorial locally, there are many directories that must be created, and many files that will need to be downloaded. In total, these files will be ~19GB in size. You will also need ~30 GB of additional free space to allow for the analysis files to be created across the tutorial steps.

**The easiest way to prepare is to download the `Directory-Structure-Data.zip` file located in the `OSF/00-Directory-Structure/` folder in the OSF files window. It will contain the complete set of directories and all the data files (except one, see below!). Please read the instructions below for more details.** 

### 1. Create Local Directories
To begin, you will need to have the following directories on your local computer:

```
bin/
│
├── 01-Starting-seqs/
├── 02-Taxon-assess/
├── 03-Rename-merge/
├── 04-Starting-materials/
├── 05-Parse-loci/
├── 06a-Cluster-blast/
│	├── Input/
│	└── Output/
├── 06b-Reference-blast/
│	├── Input/
│	└── Output/
├── 06c-Contamination-filter/
│	├── Input/
│	└── Output/
├── 07-Filter-seqs/
│	├── Input/
│	└── Output/
├── 08-Min-taxa/
├── 09-Adjust/
├── 10-Coding-tests/
│	├── Output-mtdna/
│	└── Output-nuc/
├── 11-Align/
│	├── mtDNA/
│	├── Noncoding/
│	├── Nuc/
│	└── Output/
├── 12-Acc-table/
├── 13-Relabel/
│	├── Input/
│	└── Output/
├── 14-Trim/
├── 15-Convert/
├── 16-Concatenate/
```
There are two options for getting these folders made. You can download the `Directory-Structure.zip` file or the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder and unzip it. The `Directory-Structure-Data.zip` version contains all but one of the required data files already in their correct directories (see below). The `Directory-Structure.zip` version simply contains all the empty directories. Alternatively, if you choose to not download these files you can make the directories by hand. 

Please note that the `bin/` directory is just a placeholder name, to simplify the "path" to each directory. All of the subsequent steps will reference the "paths" to these folders. These paths will be different on your local computer, so make sure to change the path as appropriate. 

In this example, the path to the `01-Starting-seqs/` folder will be called `/bin/01-Starting-seqs/`. 

### 2. Download Data

You will need to obtain several files from the OSF storage to run the analysis. The easiest way to do this is download the folders with the data already present, by getting the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder. This will have all of the files listed below, except for the starting sequences. You will need to download the starting sequences before running the tutorial!

Alternatively, if you downloaded the empty directories or made them by hand, you will need to download and move the files listed below to their correct directory.

**Starting sequences:**
The GenBank sequences downloaded in fasta format are available here in the `OSF/01-Starting-seqs/Fasta-Pieces/` folder. Because of the large size of the file (16GB), it was split into many smaller files (still totaling 16GB). If you are running through the tutorial, you will need to download these files and merge them into a single fasta file. You should call this fasta file `Iguania.fasta`. If you are using unix, you can accomplish this with a command like `cat *.fasta > Iguania.fasta` from within the folder containing the split fasta files.

When you have reconstituted the `Iguania.fasta` file, move it to your `/bin/01-Starting-seqs/` folder. 

**Taxon renaming file:** Download the `Rename_List.txt` file from `OSF/03-Rename-merge/` and move it to your `/bin/03-Rename-merge/` folder. 

**Locus search terms and taxon names list:** Download the `Locus-search-terms.txt` and `Reptile_Database_Updated_2018.txt` files from `OSF/04-Starting-materials/` and move them both to your `/bin/04-Starting-materials/` folder. 

**Files for Reference Blast:** Download all of the reference sequence files (9 total) from the `OSF/06b-Reference-blast/Input/` folder and move them to your local `/bin/06b-Reference-blast/Input/` folder. Download the `multisearch-key.txt` file from `OSF/06b-Reference-blast/` and move it to your `/bin/06b-Reference-blast/` folder.

**Files for Contamination Filtering:** Download the `Human_mtDNA.fasta` reference sequence file from the `OSF/06c-Contamination-filter/Input/` folder and move it to your local `/bin/06c-Contamination-filter/Input/` folder. Download the `multisearch-key.txt` file from `OSF/06c-Contamination-filter/` and move it to your `/bin/06c-Contamination-filter/` folder.

**Files for Sequence Filtering:** Download the three `*_File_List.txt` files from `OSF/07-Filter-seqs/` and move them to your local `/bin/07-Filter-seqs/` folder.

**Files for Translation Tests:** Download the `mtdna_fastas.txt` and `nuc_fastas.txt` files from `OSF/10-Coding-tests/` and move them to your local `/bin/10-Coding-tests/` folder. 

### 3. Check File Locations are Correct

If you downloaded the above files and moved them to the correct directories, you should find them in the following locations. The downloaded files are marked with asterisks at the beginning (***) of their names:

```
bin/
│
├── 01-Starting-seqs/
│   └── *** Iguania.fasta
├── 02-Taxon-assess/
├── 03-Rename-merge/
│   └── *** Rename_List.txt
├── 04-Starting-materials/
│   ├── *** Locus-search-terms.txt
│   └── *** Reptile_Database_Updated_2018.txt
├── 05-Parse-loci/
├── 06a-Cluster-blast/
│	├── Input/
│	└── Output/
├── 06b-Reference-blast/
│	├── *** multisearch-key.txt
│	├── Input/
│   |  ├── *** reference-sequences-vertebrates-RAG1-part1.fasta
│   |  ├── *** reference-sequences-vertebrates-RAG1-part2.fasta
│   |  ├── *** Squamate-12S-References.fasta 
│   |  ├── *** Squamate-16S-References.fasta 
│   |  ├── *** Squamate-CO1-References.fasta 
│   |  ├── *** Squamate-CYTB-References.fasta
│   |  ├── *** Squamate-ND1-References.fasta 
│   |  ├── *** Squamate-ND2-References.fasta 
│   |  └── *** Squamate-ND4-References.fasta 
│	└── Output/
├── 06c-Contamination-filter/
│	├── *** multisearch-key.txt
│	├── Input/
│   |  └── *** Human_mtDNA.fasta
│	└── Output/
├── 07-Filter-seqs/
│	├── *** mtdna_File_List.txt
│	├── *** Noncoding_File_List.txt
│	├── *** Nuc_File_List.txt
│	├── Input/
│	└── Output/
├── 08-Min-taxa/
├── 09-Adjust/
├── 10-Coding-tests/
│	├── *** mtdna_fastas.txt
│	├── *** nuc_fastas.txt
│	├── Output-mtdna/
│	└── Output-nuc/
├── 11-Align/
│	├── mtDNA/
│	├── Noncoding/
│	├── Nuc/
│	└── Output/
├── 12-Acc-table/
├── 13-Relabel/
│	├── Input/
│	└── Output/
├── 14-Trim/
├── 15-Convert/
├── 16-Concatenate/
```
If everything looks good, you are ready to proceed to the tutorial.

---------------


## **Running the Tutorial**

This tutorial will create a species-level supermatrix for Iguania. It should result in a supermatrix containing 61 loci, 1,426 species, and 13,307 total sequences. 

Just as a warning, the analysis will take ~13 hours to complete (not including user-time), due primarily to ~11.5 hours for the sequence alignment step. If you would like to run a faster tutorial, you should look at any other analysis besides this one!

The steps in this tutorial follow along with the steps listed in Supplemental Table S1 of Portik & Wiens (2020), and include:
1. Assess Taxonomy
1. Parse Loci
1. Similarity Filtering
1. Sequence Selection
1. Sequence Alignment
1. Post-Alignment

As a reminder, the paths to folders and files will likely be different on your local computer. You will need to replace the `/bin/` section of the paths used in the commands with the local path to the directories you created in your environment. 

For example, `/bin/01-Starting-seqs/` may actually be something like ` /Users/portik/Tutorial-Iguania/01-Starting-seqs/`. 


---------------

### **1. Assess Taxonomy**

#### **Taxa Assessment**

Any user-supplied taxonomy can be used in SuperCRUNCH. For this project, we want to use the current reptile taxonomy that was downloaded from the Reptile Database (file `Reptile_Database_Updated_2018.txt`). Our first step is to see if there are conflicts between the taxon names in the downloaded sequences and the taxon names from the Reptile Database. 

We will identify all name conflicts using the `Taxa_Assessment.py` module, which screens the starting sequences using the taxon names file. For this case, we are not interested in searching for subspecies names, so we turn this option off. 

Here is the command used:

```
python Taxa_Assessment.py -i /bin/01-Starting-seqs/Iguania.fasta -t /bin/04-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/02-Taxon-assess/ --no_subspecies 
```

For the unmatched and matched names, a set of accession numbers, taxon names, and a fasta file is provided as output to `/bin/02-Taxon-assess/`. 

Note that I split the output fasta files for OSF to allow them to be uploaded, but when you run this step it will output only two fasta files (`Matched_Taxa.fasta` and `Unmatched_Taxa.fasta`).

The key output is a list of the names that did not match the taxonomy of the Reptile Database (`Unmatched_Records_Species_Names.log`). I used this file to find all the taxon names that could be corrected (in general, not all names can be rescued). These names were misspelled or represented synonomies. 

I then created a file for renaming a subset of the unmatched names (`Rename_List.txt`). It contains two columns, one with an unmatched name and the other with the replacement name. This file should be located in you local `03-Rename-merge/` directory.

#### **Rename and Merge**

Given that many names can be corrected in this dataset, we next want to relabel the unmatched names in our sequence data using the `Rename_Merge.py` module.

We will relabel all sequences in the `Unmatched_Taxa.fasta` file using the corrected names provided in `Rename_List.txt`. We will then merge these updated sequence records with those contained in `Matched_Taxa.fasta`. This will produce a single, correctly labeled, output fasta file for subsequent steps.

Here is the command used:

```
python Rename_Merge.py -i /bin/02-Taxon-assess/Unmatched_Taxa.fasta -r /bin/03-Rename-merge/Rename_List.txt -m /bin/02-Taxon-assess/Matched_Taxa.fasta -o /bin/03-Rename-merge/
```

In this case, it relabeled 1,860 of 8,626,443 records in the `Unmatched_Taxa.fasta` file. This is a relatively small number, but it will increase the quality of our final dataset. 

The corrected records, plus the records from `Unmatched_Taxa.fasta`, are now present in the `Merged.fasta` file. The other records that were not relabeled (~7 million) were not included in the `Merged.fasta` file. This ended up reducing the file size of the starting sequence set from 16GB to 7.5GB for the `Merged.fasta` file. 

The `Merged.fasta` file was written to the `/bin/03-Rename-merge/` directory. 

Move the `Merged.fasta` file to the `/bin/04-Starting-materials/` folder. This will be the final set of sequences we will use to perform the main steps of SuperCRUNCH.

---------------

### **2. Parse Loci**

SuperCRUNCH will first identify loci based on the labels present in the sequence records. All sequences that match a specific gene abbreviation and/or description will be written to a fasta file created for that locus. These searches are conducted using the `Parse_Loci.py` module.

Running this step requires three files: 
- a list of genes/loci to search for
- a list of taxon names
- the fasta sequences to search. 

These files should now all be in the `/bin/04-Starting-materials/` directory, and include `Locus-search-terms.txt`, `Reptile_Database_Updated_2018.txt`, and `Merged.fasta` (from the step above). 

We will run `Parse_Loci.py` with the following command:

```
python Parse_Loci.py -i /bin/04-Starting-materials/Merged.fasta -l /bin/04-Starting-materials/Locus-search-terms.txt -t /bin/04-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/05-Parse-loci/ --no_subspecies 
```

All outputs are written to the `/bin/05-Parse-loci/` folder.

Inside there are two new folders:
- `Summary-File/` contains a text file showing the number of sequences found per locus.
- `Parsed-Fasta-Files` contains all the locus-specific fasta files for which at least two sequences were found.

The SQL database generated during this step will also be in the `/bin/05-Parse-loci/` folder.

The locus-specific fasta files contained in the `/bin/05-Parse-loci/Parsed-Fasta-Files/` folder will be used for the next step.

---------------

### **3. Similarity Filtering**

SuperCRUNCH performs similarity searches with BLASTn. It offers two options, which differ in whether the reference sequences are selected automatically or provided by the user. 

The automatic reference option is appropriate for loci containing "simple" records. We define “simple” record sets as those generally containing a single gene region with limited length variation, which results from use of the same primers (Sanger-sequencing) or probes (sequence capture) to generate sequences. 

The user-supplied reference option is appropriate for loci containing "complex" records. We define “complex” records as those containing the target region plus non-target sequence (e.g., other regions or genes). Common examples include long mtDNA fragments and whole mitogenomes, and genes that have been sequenced for different fragments that have little or no overlap. 

Sequences that have significant hits to the references will be trimmed if necessary. This means that only the target regions of the input sequences are included in the output fasta file. 

#### **Run Cluster_Blast_Extract for nuclear loci ("simple" records)**

Similarity searches using automatically selected references can be performed with the `Cluster_Blast_Extract.py` module. We will use this module to filter all the nuclear loci (except for one). 

We will need to copy or move the fasta files for the nuclear loci from the `/bin/05-Parse-loci/Parsed-Fasta-Files/` folder to the `/bin/06a-Cluster-blast/Input/` folder. 

Copy or move ALL the files in the `/bin/05-Parse-loci/Parsed-Fasta-Files/` folder EXCEPT for the following 8 files:

`12S.fasta`, `16S.fasta`, `CO1.fasta`, `CYTB.fasta`, `ND1.fasta`, `ND2.fasta`, `ND4.fasta`, `RAG1.fasta`. 

This should result in 58 files present in the `/bin/06a-Cluster-blast/Input/` folder.

We are now ready to run the `Cluster_Blast_Extract.py` module for these 58 nuclear loci.

We will run the module with the following command:

```
python Cluster_Blast_Extract.py -i /bin/06a-Cluster-blast/Input/ -o /bin/06a-Cluster-blast/Output/ -b dc-megablast --max_hits 200 -m span --threads 4
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Four new directories have been written to the `/bin/06a-Cluster-blast/Output/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/06a-Cluster-blast/Output/Filtered-Fasta-Files/` directory. We will use these fasta files for the sequence selection step. 


#### **Run Reference_Blast_Extract for loci with "complex" records**

Similarity searches using user-supplied references can be performed with the `Reference_Blast_Extract.py` module. We will use this module to filter the mtDNA genes, plus the nuclear gene RAG-1. In reptiles, RAG-1 has historically been sequenced for two non-overlapping fragments, so we will search for both fragments separately.

We will need to copy or move the fasta files for these loci from the `/bin/05-Parse-loci/Parsed-Fasta-Files/` folder to the `/bin/06b-Reference-blast/Input/` folder. 

Copy or move the following 8 files located in `/bin/05-Parse-loci/Parsed-Fasta-Files/`:

`12S.fasta`, `16S.fasta`, `CO1.fasta`, `CYTB.fasta`, `ND1.fasta`, `ND2.fasta`, `ND4.fasta`, `RAG1.fasta`. 

Because RAG1 will need to be screened twice, we will duplicate the file. Make a copy of `RAG1.fasta`, and rename the two RAG1 files to be: `RAG1p1.fasta` and `RAG1p2.fasta`.

This should result in 9 files present in the `/bin/06b-Reference-blast/Input/` folder.

In addition to the fasta files, there should also be 9 corresponding reference sequence files (see Preparing for the Tutorial) in the `/bin/06b-Reference-blast/Input/` folder.

We will take advantage of the multisearch option of `Reference_Blast_Extract.py` to automate runs. This will pair each fasta file with a reference file and run them sequentially. To use this option we will need the `multisearch-key.txt` file, which should be located in the `/bin/06b-Reference-blast/` folder.

We are now ready to run the `Reference_Blast_Extract.py` module for these 9 loci. 

We will run the module with the following command:

```
python Reference_Blast_Extract.py -i /bin/06b-Reference-blast/Input/ -o /bin/06b-Reference-blast/Output/ --multisearch /bin/06b-Reference-blast/multisearch-key.txt -b dc-megablast --max_hits 200 -m span --threads 4
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Three new directories have been written to the `/bin/06b-Reference-blast/Output/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/06b-Reference-blast/Output/Filtered-Fasta-Files/` directory. 

The filtered RAG1 fasta files are ready for the sequence selection step.

Before sequence selection, we will run the mtDNA loci through one additional filter, which is described below.

#### **Run Contamination_Filter for the filtered mtDNA loci**

The contamination filter can be used to remove all sequences scoring >95% identity for at least 100 continuous base pairs to a particular reference sequence (the source of contamination). This is achieved using the `Contamination_Filter.py` module. 

For our use case, we want to screen the mtDNA loci (which have already passed through similarity filtering) against human mtDNA sequences.

To prepare, you'll need to copy the filtered mtDNA fasta files from the `/bin/06b-Reference-blast/Output/Filtered-Fasta-Files/` to the input bin for this step, which is `/bin/06c-Contamination-filter/Input`. 

There are 7 mtDNA files that need to be copied: `12S_extracted.fasta`, `16S_extracted.fasta`, `CO1_extracted.fasta`, `CYTB_extracted.fasta`, `ND1_extracted.fasta`, `ND2_extracted.fasta`, `ND4_extracted.fasta`.

In addition to the fasta files, there should also be one contamination reference sequence file called `Human_mtDNA.fasta` (see Preparing for the Tutorial) in the `/bin/06c-Contamination-filter/Input/` folder.

We will take advantage of the multisearch option of `Contamination_Filter.py` to automate runs. This will pair each fasta file to the human mtDNA reference file and run them sequentially. To use this option we will need the `multisearch-key.txt` file, which should be located in the `/bin/06c-Contamination-filter/` folder.

We will run the contamination filter using the following command:

```
python Contamination_Filter.py -i /bin/06c-Contamination-filter/Input/ -o /bin/06c-Contamination-filter/Output --multisearch /bin/06c-Contamination-filter/multisearch-key.txt -b megablast 
```

Three new directories have been written to the `/bin/06c-Contamination-filter/Output/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/06c-Contamination-filter/Output/Filtered-Fasta-Files/` directory. 

Looking in the other output directories, we can see that the contamination filter removed one 12S sequence and one 16S sequence. An online BLAST using NCBI confirms that these are indeed human sequences.

The filtered mtDNA fasta files are ready for the sequence selection step.

---------------

### **4. Sequence Selection**

For supermatrices, a single sequence is typically used to represent each species for each gene. If multiple sequences for a given gene exist for a given species (e.g., because multiple individuals were sampled), then an objective strategy must be used for sequence selection. We can use the `Filter_Seqs_and_Species.py` module to accomplish this. 

We are going to select sequences slightly differently for nuclear loci, mtDNA coding loci, and mtDNA non-coding loci. 

To prepare for all the subsequent steps, you will need to move or copy ALL of the filtered fasta files from the previous step to the `/bin/07-Filter-seqs/Input/` folder.

There should be 67 filtered fasta files in total. This includes:

- the 58 filtered nuclear loci files from: `/bin/06a-Cluster-blast/Output/Filtered-Fasta-Files/`
- the 2 filtered RAG1 files from: `/bin/06b-Reference-blast/Output/Filtered-Fasta-Files/`
- the 7 contamination-filtered mtDNA loci files from: `/bin/06c-Contamination-filter/Output/Filtered-Fasta-Files/`

We will use the file inclusion option of `Filter_Seqs_and_Species.py` to run analyses for different sets of fasta files. To take advantage of this, we need to ensure that three files are present in the `/bin/07-Filter-seqs/` folder: `mtdna_File_List.txt`, `Noncoding_File_List.txt`, `Nuc_File_List.txt`.

These three files should already be present in this folder (see Preparing for the Tutorial). 

We are now ready to run each of the three sequence selection steps.

#### **Run Filter_Seqs_and_Species for protein-coding nuclear loci**

We will select sequences for the nuclear protein-coding loci based on length (longest) and whether the sequences have a valid reading frame (they can be translated, no errors/stop codons). We will specify to translate sequences using standard code, and set 200 bp as the minimum length required to keep a sequence.

To run this step, we need the filtered fasta files, a taxon names list file (`Reptile_Database_Updated_2018.txt`), and the file specifying which fasta files to include (`Nuc_File_List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/07-Filter-seqs/Input/ -s oneseq -f translate --table standard -m 200 -t /bin/04-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/07-Filter-seqs/Output/ --onlyinclude /bin/07-Filter-seqs/Nuc_File_List.txt --no_subspecies --quiet 
```

The outputs are written to three folders contained in `/bin/07-Filter-seqs/Output/Results/`. The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/07-Filter-seqs/Output/Results/Filtered-Fasta-Files/`.

#### **Run Filter_Seqs_and_Species for protein-coding mtDNA loci**

We will select sequences for the mtDNA protein-coding loci based on length (longest) and whether the sequences have a valid reading frame (they can be translated, no errors/stop codons). We will specify to translate sequences using vertebrate mtDNA code, and set 200 bp as the minimum length required to keep a sequence.

To run this step, we need the filtered fasta files, a taxon names list file (`Reptile_Database_Updated_2018.txt`), and the file specifying which fasta files to include (`mtdna_File_List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/07-Filter-seqs/Input/ -s oneseq -f translate --table vertmtdna -m 200 -t /bin/04-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/07-Filter-seqs/Output/ --onlyinclude /bin/07-Filter-seqs/mtdna_File_List.txt --no_subspecies --quiet 
```

The outputs are written to three folders contained in `/bin/07-Filter-seqs/Output/Results/`. The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/07-Filter-seqs/Output/Results/Filtered-Fasta-Files/`.


#### **Run Filter_Seqs_and_Species for non-protein-coding mtDNA loci**

We will select sequences for the mtDNA non-protein-coding loci based only on length (longest). We will set 200 bp as the minimum length required to keep a sequence.

To run this step, we need the filtered fasta files, a taxon names list file (`Reptile_Database_Updated_2018.txt`), and the file specifying which fasta files to include (`Noncoding_File_List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/07-Filter-seqs/Input/ -s oneseq -f length -m 200 -t /bin/04-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/07-Filter-seqs/Output/ --onlyinclude /bin/07-Filter-seqs/Noncoding_File_List.txt --no_subspecies --quiet 
```

The outputs are written to three folders contained in `/bin/07-Filter-seqs/Output/Results/`. The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/07-Filter-seqs/Output/Results/Filtered-Fasta-Files/`.


#### **Filter Loci by a minimum sequence requirement**

Before moving on to alignment, we will remove all loci that have fewer than 30 sequences. For this analysis, the sequences now represent species (we just selected one sequence per species). So, we want to ensure that at least 30 species are present for each locus.

To apply this filter, we can use the `Fasta_Filter_by_Min_Seqs.py` module.

All of the fasta files resulting from sequence selection are in the `/bin/07-Filter-seqs/Output/Results/Filtered-Fasta-Files/` folder.

We will run `Fasta_Filter_by_Min_Seqs.py` using the following command:

```
python Fasta_Filter_by_Min_Seqs.py -i /bin/07-Filter-seqs/Output/Results/Filtered-Fasta-Files/ --min_seqs 30 -o /bin/08-Min-taxa/
```

The output on screen tells us that 61 fasta files passed the minimum sequence filter, and 6 fasta files failed the minimum sequence filter.

The fasta files that passed the filter will now be in the `/bin/08-Min-taxa/Filtered-Fasta-Files/` directory.


---------------

### **5. Sequence Alignment**

The overall strategy is to perform translation alignments for the nuclear protein-coding loci and the mtDNA protein-coding loci. This will be done using the alignment program MACSE. We will use Clustal-O to align non-protein-coding mtDNA loci.  

In preparation for alignment, we need to complete three pre-alignment steps. 

First, we will adjust the direction of all sequences within each fasta file. This will ensure they are all oriented in the same direction, and prevent catastrophes in sequence alignment.

Second, we will attempt to translate all nuclear protein-coding loci into their correct reading frames, and adjust them to the first codon position. This step will produce all the input files necessary for performing translation alignments with the alignment program MACSE.

Third, we will attempt to translate all mtDNA protein-coding loci into their correct reading frames, and adjust them to the first codon position. This step will produce all the input files necessary for performing translation alignments with the alignment program MACSE.

To prepare for the pre-alignment steps, please ensure that `mtdna_fastas.txt` and `nuc_fastas.txt` are both in the `/bin/10-Coding-tests/` directory.

#### **Pre-Alignment: Adjust sequence directions**

We will ensure the sequences in each fasta file are correctly oriented by using the `Adjust_Direction.py` module. All of the fasta files should be in the same output directory after the minimum sequence filter, so we will just run this step using that location. 

Here is the command we will use:

```
python Adjust_Direction.py -i /bin/08-Min-taxa/Filtered-Fasta-Files/ -o /bin/09-Adjust/ --threads 8
```

You should adjust the number of threads accordingly.

There should now be two output directories and a summary file in the `/bin/09-Adjust/` folder. The adjusted fasta files are the main output we are interested in, and they are located in the `/bin/09-Adjust/Adjusted-Fasta-Files/` folder.

#### **Pre-Alignment: Translate nuclear protein-coding loci**

We will now translate and adjust reading frames for the nuclear protein-coding loci. We will do this using the `Coding_Translation_Tests.py` module, using standard code for translation. We will use the file inclusion feature to only select a subset of the fasta files from the `/bin/09-Adjust/Adjusted-Fasta-Files/` folder.

We will run this step with the following command:

```
python Coding_Translation_Tests.py -i /bin/09-Adjust/Adjusted-Fasta-Files/ -o /bin/10-Coding-tests/Output-nuc/ --table standard --onlyinclude /bin/10-Coding-tests/nuc_fastas.txt --quiet
```

Three output directories are created in the `/bin/10-Coding-tests/Output-nuc/`. For MACSE, we will use the files in the `Translation-Failed-Seqs/` and `Translation-Passed-Seqs/` folders.


#### **Pre-Alignment: Translate mtDNA protein-coding loci**

We will now translate and adjust reading frames for the mtDNA protein-coding loci. We will do this using the `Coding_Translation_Tests.py` module, using vertebrate mitochondrial code for translation. We will use the file inclusion feature to only select a subset of the fasta files from the `/bin/09-Adjust/Adjusted-Fasta-Files/` folder.

We will run this step with the following command:

```
python Coding_Translation_Tests.py -i /bin/09-Adjust/Adjusted-Fasta-Files/ -o /bin/10-Coding-tests/Output-mtdna/ --table vertmtdna --onlyinclude /bin/10-Coding-tests/mtdna_fastas.txt --quiet
```

Three output directories are created in the `/bin/10-Coding-tests/Output-mtdna/`. For MACSE, we will use the files in the `Translation-Failed-Seqs/` and `Translation-Passed-Seqs/` folders.

#### **Alignment: Nuclear protein-coding loci**

We will now use the translation alignment method in MACSE to align the protein-coding nuclear loci. This step will take at least ~6 hours to complete. 

To prepare, copy or move ALL of the fasta files present in the `/bin/10-Coding-tests/Output-nuc/Translation-Failed-Seqs/` and `/bin/10-Coding-tests/Output-nuc/Translation-Passed-Seqs/` folders to the `/bin/11-Align/Nuc/` folder. 

We will use the algorithm in MACSE to correct reading frame errors in the sequences that failed translation, using all of the sequences that passed translation, for each locus. Thus, we will use the pass_fail option. 

You will need to have the MACSE (v2) jar file downloaded on your computer, and must specify the full path to this file. Here, it will just say `/bin/macse/macse_v2.03.jar`. 

We will perform these alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/11-Align/Nuc/ -o /bin/11-Align/Output/ -a macse --mpath /bin/macse/macse_v2.03.jar --table standard --mem 10 --pass_fail
```

Note that the memory is set to 10GB here, please adjust accordingly.

The final alignments will be written to `/bin/11-Align/Output/Alignments-MACSE/Cleaned-Alignments/`.


#### **Alignment: mtDNA protein-coding loci**

We will now use the translation alignment method in MACSE to align the protein-coding mtDNA loci. This step will take at least ~6 hours to complete. 

To prepare, copy or move ALL of the fasta files present in the `/bin/10-Coding-tests/Output-mtdna/Translation-Failed-Seqs/` and `/bin/10-Coding-tests/Output-mtdna/Translation-Passed-Seqs/` folders to the `/bin/11-Align/mtDNA/` folder. 

We will use the algorithm in MACSE to correct reading frame errors in the sequences that failed translation, using all of the sequences that passed translation, for each locus. Thus, we will use the pass_fail option. 

You will need to have the MACSE (v2) jar file downloaded on your computer, and must specify the full path to this file. Here, it will just say `/bin/macse/macse_v2.03.jar`. 

We will perform these alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/11-Align/mtDNA/ -o /bin/11-Align/Output/ -a macse --mpath /bin/macse/macse_v2.03.jar --table vertmtdna --mem 10 --pass_fail
```

Note that the memory is set to 10GB here, please adjust accordingly.

**Warning:** Sometimes MACSE seemed buggy with respect to certain loci, and the alignments failed. This seemed to happen with CO1 and CYTB. If this happens for you, try using Clustal-O instead.

The final alignments will be written to `/bin/11-Align/Output/Alignments-MACSE/Cleaned-Alignments/`.

#### **Alignment: mtDNA non-protein-coding loci**

We will use Clustal-O to align the remaining loci, which include 12S and 16S. This step will take ~30 min.

To prepare, please move or copy the `12S_extracted_contfiltered_oneseq_Adjusted.fasta` and `16S_extracted_contfiltered_oneseq_Adjusted.fasta` files from the `/bin/09-Adjust/Adjusted-Fasta-Files/` folder to the `/bin/11-Align/Noncoding/` folder. 

We will perform these alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/11-Align/Noncoding/ -o /bin/11-Align/Output/ -a clustalo --accurate --threads 4
```

The final alignments will be written to `/bin/11-Align/Output/Alignments-CLUSTALO/`.


---------------

### **6. Post-Alignment**

We will now perform several post-alignment tasks:

- We will create a table of GenBank accession numbers for all the sequences included in our alignments.
- We will relabel the sequence records with relevant taxon names, making them compatible with downstream programs.
- We will trim alignments.
- We will convert the fasta alignments into nexus and phylip format.
- We will concatenate the alignments to make a final supermatrix that is compatible with RAxML.

#### **Make a table of GenBank accession numbers**

Organizing GenBank accession numbers can be a serious undertaking, but this can be automated for you by using the `Make_Acc_Table.py` module. 

This step can be performed any time after sequence selection, but is best to use after alignment. For this example, we will use the files generated right after the minimum sequence requirement filter (before alignment).

We can generate this handy table using the following command:

```
python Make_Acc_Table.py -i /bin/08-Min-taxa/Filtered-Fasta-Files/ -o /bin/12-Acc-table/
```

The outputs include the table itself (`GenBank_Accession_Table.txt`) and a summary of the taxa included and the number of sequences/loci for each taxon (`Taxon_Summary.txt`).


#### **Relabel sequences in alignments**

Currently, the alignments contain sequences that are labeled using a modified version of the GenBank record line. This cannot be used for most downstream programs. We will relabel these sequences to include only the taxon label.

To run this step, we will need to move all of the alignment files from the `/bin/11-Align/Output/Alignments-CLUSTALO/` and `/bin/11-Align/Output/Alignments-MACSE/Cleaned-Alignments` folders to the `/bin/13-Relabel/Input/` folder. There should be 61 alignment files in total.

We will relabel sequences with the `Fasta_Relabel_Seqs.py` module using the following command:

```
python Fasta_Relabel_Seqs.py -i /bin/13-Relabel/Input/ -o /bin/13-Relabel/Output/ -r species
```

The fasta alignments with relabeled sequences will be written to `/bin/13-Relabel/Output/Relabeled-by-species/`.

#### **Trim alignments**

We will perform light trimming to remove very poorly aligned regions from all the alignments. We will use the `Trim_Alignments_Trimal.py` module to accomplish this. We will use the relabeled fasta alignments found in `/bin/13-Relabel/Output/Relabeled-by-species/`.


We will trim alignments using the following command:

```
python Trim_Alignments_Trimal.py -i /bin/13-Relabel/Output/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/14-Trim/ -f fasta -a gt --gt 0.1
```

This will run the gap-threshold algorithm of trimal, using a threshold value of 0.1. 

The trimmed alignments will be written to the `/bin/14-Trim/Trimmed-fasta-Files/` directory.

#### **Convert alignments into nexus and phylip format**

We can quickly convert any relabeled fasta alignment to nexus and phylip format using the `Fasta_Convert.py` module. We will perform this format conversion for all of the trimmed alignments found in `/bin/14-Trim/Trimmed-fasta-Files/`.

We will convert formats using the following command:

```
python Fasta_Convert.py -i /bin/14-Trim/Trimmed-fasta-Files/ -o /bin/15-Convert/
```

The nexus-formatted files will be written to `/bin/15-Convert/Nexus-Files/` and the phylip-formatted files will be written to `/bin/15-Convert/Phylip-Files/`.

#### **Concatenate alignments to create the final supermatrix**

The moment is upon us! It is time to create the final supermatrix. We will use the `Concatenation.py` module to perform this task. We will use the trimmed fasta files for this step. 

We will specify the input format is fasta, the output format is phylip (compatible with RAxML), and that missing characters should be coded with a `-` symbol.

We will perform this task using the following command:

```
python Concatenation.py -i /bin/14-Trim/Trimmed-fasta-Files/ -o /bin/16-Concatenate/ --informat fasta --outformat phylip -s dash
```

And here is the output on screen when we run it:

```
Found 61 fasta files to concatenate.

Found 1,426 unique taxa across alignment files.

Gathering sequences for all taxa (this could take some time)...
	Done.

Writing partitions file.

Writing concatenated alignment.
	Total alignment length = 53,319 bp.
	Total number of sequences included = 13,307.

Finished. Total elapsed time: 0:00:00.630740 (H:M:S)
```

A concatenated alignment in phylip format (`Concatenated_Alignment.phy`), a data partitions file (`Data_Partitions.txt`), and final count of the number of taxa and their number of corresponding loci (`Taxa_Loci_Count.log`) will be written to `/bin/16-Concatenate/`. 

And that's it! We've just created the most comprehensive and high-quality supermatrix ever assembled for iguanian lizards. 

You can use the concatenated alignment to create a phylogeny using RAxML, like we did in our study. Our bootstrapped tree is available in the `OSF/17-RAxML/` folder, in case you'd like to download it.


