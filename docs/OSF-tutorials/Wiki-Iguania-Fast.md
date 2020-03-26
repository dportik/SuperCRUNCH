## Tutorial for Analysis: Iguania-Fast

This wiki is intended to allow you to replicate the Iguania-Fast analysis reported in the Supplemental Materials. All of the input and output files of the completed analysis are provided in the OSF data folders.  

The wiki is presented as a tutorial. However, you can always scroll through the sections of this wiki to see how each module was called, and which files were used and produced at each step. One of the main goals is to provide you with all the input files required for running the different modules of SuperCRUNCH, and to demonstrate the use of various options along the way. This page is therefore intended to function as a tutorial and a reference resource.

**Please note that this analysis requires ~15GB of space and will take at least ~1.5 hours to complete.**

---------------


## **Preparing for the Tutorial**

Before you can run this tutorial locally, there are many directories that must be created, and many files that will need to be downloaded. In total, these files will be ~13GB in size. You will also need ~2GB of additional free space to allow for the analysis files to be created across the tutorial steps.

**The easiest way to prepare is to download the `Directory-Structure-Data.zip` file located in the `OSF/00-Directory-Structure/` folder in the OSF files window. It will contain the complete set of directories and all the data files (except one, see below!). Please read the instructions below for more details.** 

### 1. Create Local Directories
To begin, you will need to have the following directories on your local computer:

```
bin/
│
├── 01-Starting-materials/
├── 02-Parse-loci/
├── 03a-Cluster-blast/
│	├── Input/
│	└── Output/
├── 03b-Reference-blast/
│	├── Input/
│	└── Output/
├── 04-Filter-seqs/
│	├── Input/
│	└── Output/
├── 05-Min-taxa/
├── 06-Make-acc-table/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Trim/
├── 11-Concatenate/
```
There are two options for getting these folders made. You can download the `Directory-Structure.zip` file or the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder and unzip it. The `Directory-Structure-Data.zip` version contains all but one of the required data files already in their correct directories (see below). The `Directory-Structure.zip` version simply contains all the empty directories. Alternatively, if you choose to not download these files you can make the directories by hand. 

Please note that the `bin/` directory is just a placeholder name, to simplify the "path" to each directory. All of the subsequent steps will reference the "paths" to these folders. These paths will be different on your local computer, so make sure to change the path as appropriate. 

In this example, the path to the `01-Starting-materials/` folder will be called `/bin/01-Starting-materials/`. 

### 2. Download Data

You will need to obtain several files from the OSF storage to run the analysis. The easiest way to do this is download the folders with the data already present, by getting the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder. This will have all of the files listed below, except for the starting sequences. You will need to download the starting sequences before running the tutorial!

Alternatively, if you downloaded the empty directories or made them by hand, you will need to download and move the files listed below from their OSF directories and move them to their correct local directory.

**Starting sequences:**
The GenBank sequences downloaded in fasta format are available here in the `OSF/01-Starting-materials/Starting-seqs/Fasta-Pieces/` folder. Because of the large size of the file (16GB), it was split into many smaller files (still totaling 16GB). If you are running through the tutorial, you will need to download these files and merge them into a single fasta file. You should call this fasta file `Iguania.fasta`. If you are using unix, you can accomplish this with a command like `cat *.fasta > Iguania.fasta` from within the folder containing the split fasta files.

When you have reconstituted the `Iguania.fasta` file, move it to your local `/bin/01-Starting-materials/` folder. 

**Locus search terms and taxon names list:** Download the `Locus-search-terms.txt` and `Reptile_Database_Updated_2018.txt` files from `OSF/01-Starting-materials/` and move them both to your `/bin/01-Starting-materials/` folder. 

**Files for Reference Blast:** Download all of the reference sequence files (9 total) from the `OSF/03b-Reference-blast/Input/` folder and move them to your local `/bin/03b-Reference-blast/Input/` folder. Download the `multisearch-key.txt` file from `OSF/03b-Reference-blast/` and move it to your `/bin/03b-Reference-blast/` folder.

### 3. Check File Locations are Correct

If you downloaded the above files and moved them to the correct directories, you should find them in the following locations. The downloaded files are marked with asterisks at the beginning (***) of their names:

```
bin/
│
├── 01-Starting-materials/
│   ├── *** Iguania.fasta
│   ├── *** Locus-search-terms.txt
│   └── *** Reptile_Database_Updated_2018.txt
├── 02-Parse-loci/
├── 03a-Cluster-blast/
│	├── Input/
│	└── Output/
├── 03b-Reference-blast/
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
├── 04-Filter-seqs/
│	├── Input/
│	└── Output/
├── 05-Min-taxa/
├── 06-Make-acc-table/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Trim/
├── 11-Concatenate/
```

If everything looks good, you are ready to proceed to the tutorial.

---------------


## **Running the Tutorial**

This tutorial will create a species-level supermatrix for Iguania. It should result in a supermatrix containing 60 loci, 1,399 species, and 12,978 total sequences. 

This analysis will take roughly ~1.5 hours to run (not including user-time), and is a good example of how to create a large supermatrix quickly.

The steps in this tutorial follow along with the steps listed in Supplemental Table S2 of Portik & Wiens (2020), and include:
1. Parse Loci
1. Similarity Filtering
1. Sequence Selection
1. Sequence Alignment
1. Post-Alignment

As a reminder, the paths to folders and files will likely be different on your local computer. You will need to replace the `/bin/` section of the paths used in the commands with the local path to the directories you created in your environment. 

For example, `/bin/01-Starting-materials/` may actually be something like ` /Users/portik/Tutorial-Iguania/01-Starting-materials/`. 


---------------

### **1. Parse Loci**

SuperCRUNCH will first identify loci based on the labels present in the sequence records. All sequences that match a specific gene abbreviation and/or description will be written to a fasta file created for that locus. These searches are conducted using the `Parse_Loci.py` module.

Running this step requires three files: 
- a list of genes/loci to search for
- a list of taxon names
- the fasta sequences to search. 

These files should now all be in the `/bin/01-Starting-materials/` directory, and include `Locus-search-terms.txt`, `Reptile_Database_Updated_2018.txt`, and `Iguania.fasta`. 

We will run `Parse_Loci.py` with the following command:

```
python Parse_Loci.py -i /bin/01-Starting-Materials/Iguania.fasta -l /bin/01-Starting-materials/Locus-search-terms.txt -t /bin/01-Starting-materials/Reptile_Database_Updated_2018.txt -o /bin/02-Parse-loci/ --no_subspecies 
```

All outputs are written to the `/bin/02-Parse-loci/` folder.

Inside there are two new folders:
- `Summary-File/` contains a text file showing the number of sequences found per locus.
- `Parsed-Fasta-Files` contains all the locus-specific fasta files for which at least two sequences were found.

The SQL database generated during this step will also be in the `/bin/02-Parse-loci/` folder.

The locus-specific fasta files contained in the `/bin/02-Parse-loci/Parsed-Fasta-Files/` folder will be used for the next step.

---------------

### **2. Similarity Filtering**

SuperCRUNCH performs similarity searches with BLASTn. It offers two options, which differ in whether the reference sequences are selected automatically or provided by the user. 

The automatic reference option is appropriate for loci containing "simple" records. We define “simple” record sets as those generally containing a single gene region with limited length variation, which results from use of the same primers (Sanger-sequencing) or probes (sequence capture) to generate sequences. 

The user-supplied reference option is appropriate for loci containing "complex" records. We define “complex” records as those containing the target region plus non-target sequence (e.g., other regions or genes). Common examples include long mtDNA fragments and whole mitogenomes, and genes that have been sequenced for different fragments that have little or no overlap. 

Sequences that have significant hits to the references will be trimmed if necessary. This means that only the target regions of the input sequences are included in the output fasta file. 

#### **Run Cluster_Blast_Extract for loci with "simple" records**

Similarity searches using automatically selected references can be performed with the `Cluster_Blast_Extract.py` module. We will use this module to filter all the nuclear loci (except for one). 

We will need to copy or move the fasta files for the nuclear loci from the `/bin/02-Parse-loci/Parsed-Fasta-Files/` folder to the `/bin/03a-Cluster-blast/Input/` folder. 

Copy or move ALL the files in the `/bin/02-Parse-loci/Parsed-Fasta-Files/` folder EXCEPT for the following 8 files:

`12S.fasta`, `16S.fasta`, `CO1.fasta`, `CYTB.fasta`, `ND1.fasta`, `ND2.fasta`, `ND4.fasta`, `RAG1.fasta`. 

This should result in 58 files present in the `/bin/03a-Cluster-blast/Input/` folder.

We are now ready to run the `Cluster_Blast_Extract.py` module for these 58 nuclear loci.

We will run the module with the following command:

```
python Cluster_Blast_Extract.py -i /bin/03a-Cluster-blast/Input/ -o /bin/03a-Cluster-blast/Output/ -b dc-megablast --max_hits 200 -m span --threads 4
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Four new directories have been written to the `/bin/03a-Cluster-blast/Output/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/03a-Cluster-blast/Output/Filtered-Fasta-Files/` directory. We will use these fasta files for the sequence selection step. 


#### **Run Reference_Blast_Extract for loci with "complex" records**

Similarity searches using user-supplied references can be performed with the `Reference_Blast_Extract.py` module. We will use this module to filter the mtDNA genes, plus the nuclear gene RAG-1. In reptiles, RAG-1 has historically been sequenced for two non-overlapping fragments, so we will search for both fragments separately.

We will need to copy or move the fasta files for these loci from the `/bin/02-Parse-loci/Parsed-Fasta-Files/` folder to the `/bin/03b-Reference-blast/Input/` folder. 

Copy or move the following 8 files located in `/bin/02-Parse-loci/Parsed-Fasta-Files/`:

`12S.fasta`, `16S.fasta`, `CO1.fasta`, `CYTB.fasta`, `ND1.fasta`, `ND2.fasta`, `ND4.fasta`, `RAG1.fasta`. 

Because RAG1 will need to be screened twice, we will duplicate the file. Make a copy of `RAG1.fasta`, and rename the two RAG1 files to be: `RAG1p1.fasta` and `RAG1p2.fasta`.

This should result in 9 files present in the `/bin/03b-Reference-blast/Input/` folder.

In addition to the fasta files, there should also be 9 corresponding reference sequence files (see Preparing for the Tutorial) in the `/bin/03b-Reference-blast/Input/` folder.

We will take advantage of the multisearch option of `Reference_Blast_Extract.py` to automate runs. This will pair each fasta file with a reference file and run them sequentially. To use this option we will need the `multisearch-key.txt` file, which should be located in the `/bin/03b-Reference-blast/` folder.

We are now ready to run the `Reference_Blast_Extract.py` module for these 9 loci. 

We will run the module with the following command:

```
python Reference_Blast_Extract.py -i /bin/03b-Reference-blast/Input/ -o /bin/03b-Reference-blast/Output/ --multisearch /bin/03b-Reference-blast/multisearch-key.txt -b dc-megablast -m span --threads 4 --max_hits 200
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Three new directories have been written to the `/bin/03b-Reference-blast/Output/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/03b-Reference-blast/Output/Filtered-Fasta-Files/` directory. 

These filtered fasta files are ready for the sequence selection step.

---------------

### **3. Sequence Selection**

Here we select a single sequence per species per locus, enforce a minimum number of taxa per locus, and make a table of GenBank accession numbers.

#### **Select Representative Sequences**

For supermatrices, a single sequence is typically used to represent each species for each gene. If multiple sequences for a given gene exist for a given species (e.g., because multiple individuals were sampled), then an objective strategy must be used for sequence selection. We can use the `Filter_Seqs_and_Species.py` module to accomplish this. 

We will select sequences for all loci based only on length (longest). We will also set 200 bp as the minimum length required to keep a sequence.


To prepare for this step, you will need to move or copy ALL of the filtered fasta files from the previous step to the `/bin/04-Filter-seqs/Input/` folder.

There should be 67 filtered fasta files in total. This includes:

- the 58 filtered fasta files from: `/bin/03a-Cluster-blast/Output/Filtered-Fasta-Files/`
- the 9 filtered fasta files from: `/bin/03b-Reference-blast/Output/Filtered-Fasta-Files/`

To run `Filter_Seqs_and_Species.py`, we need the filtered fasta files and a taxon names list file (`Reptile_Database_Updated_2018.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/04-Filter-seqs/Input/ -o /bin/04-Filter-seqs/Output/ -s oneseq -f length -m 200 -t /bin/01-Starting-materials/Reptile_Database_Updated_2018.txt --no_subspecies --quiet
```

The outputs are written to three folders contained in `/bin/04-Filter-seqs/Output/Results/`. 

The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/04-Filter-seqs/Output/Results/Filtered-Fasta-Files/`. We will use these for the next step.

#### **Filter Loci by a minimum sequence requirement**

Before moving on to alignment, we will remove all loci that have fewer than 30 sequences. For this analysis, the sequences now represent species (we just selected one sequence per species). So, we want to ensure that at least 30 species are present for each locus.

To apply this filter, we can use the `Fasta_Filter_by_Min_Seqs.py` module.

All of the fasta files resulting from sequence selection are in the `/bin/04-Filter-seqs/Output/Results/Filtered-Fasta-Files/` folder.

We will run `Fasta_Filter_by_Min_Seqs.py` using the following command:

```
python Fasta_Filter_by_Min_Seqs.py -i /bin/04-Filter-seqs/Output/Results/Filtered-Fasta-Files --min_seqs 30 -o /bin/05-Min-taxa 
```

The output on screen tells us that 60 fasta files passed the minimum sequence filter, and 7 fasta files failed the minimum sequence filter.

The fasta files that passed the filter will now be in the `/bin/05-Min-taxa/Filtered-Fasta-Files` directory.

#### **Make a table of GenBank accession numbers**

Organizing GenBank accession numbers can be a serious undertaking, but this can be automated for you by using the `Make_Acc_Table.py` module. 

This step can be performed any time after sequence selection, and here we perform it after the minimum taxon filter.

We can generate this handy table using the following command:

```
python Make_Acc_Table.py -i /bin/05-Min-taxa/Filtered-Fasta-Files/ -o /bin/06-Make-acc-table/
```

The outputs include the table itself (`GenBank_Accession_Table.txt`) and a summary of the taxa included and the number of sequences/loci for each taxon (`Taxon_Summary.txt`).

---------------

### **4. Sequence Alignment**

In preparation for alignment, we need to complete one pre-alignment step. We will adjust the direction of all sequences within each fasta file. This will ensure they are all oriented in the same direction, and prevent catastrophes in sequence alignment.


#### **Pre-Alignment: Adjust sequence directions**

We will ensure the sequences in each fasta file are correctly oriented by using the `Adjust_Direction.py` module. All of the target fasta files should be in the same output directory after the minimum sequence filter, so we will just run this step using that location. 

Here is the command we will use:

```
python Adjust_Direction.py -i /bin/05-Min-taxa/Filtered-Fasta-Files/ -o /bin/07-Adjust/ --threads 8
```

You should adjust the number of threads accordingly.

There should now be two output directories and a summary file in the `/bin/07-Adjust/` folder. The adjusted fasta files are the main output we are interested in, and they are located in the `/bin/07-Adjust/Adjusted-Fasta-Files/` folder. We will use these for alignment.

#### **Multiple Sequence Alignment**

We will use MAFFT (--auto) to align all loci. This is by far the fastest alignment option, and should take ~5 minutes total.

We will perform sequence alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/07-Adjust/Adjusted-Fasta-Files/ -o /bin/08-Align/ -a mafft --threads 8
```

You should adjust the number of threads accordingly.

The final alignments will be written to `/bin/08-Align/Output/Alignments-MAFFT/`.

---------------

### **5. Post-Alignment**

We will now perform several post-alignment tasks:

- We will relabel the sequence records with relevant taxon names, making them compatible with downstream programs.
- We will trim alignments.
- We will concatenate the alignments to make a final supermatrix.

#### **Relabel sequences in alignments**

Currently, the alignments contain sequences that are labeled using a modified version of the GenBank record line. This cannot be used for most downstream programs. We will relabel these sequences to include only the taxon label.

All of the target fasta files are in the `/bin/08-Align/Alignments-MAFFT/` folder.

We will relabel sequences with the `Fasta_Relabel_Seqs.py` module using the following command:

```
python Fasta_Relabel_Seqs.py -i /bin/08-Align/Alignments-MAFFT/ -o /bin/09-Relabel/ -r species 
```

The fasta alignments with relabeled sequences will be written to `/bin/09-Relabel/Output/Relabeled-by-species/`.

#### **Trim alignments**

We will perform light trimming to remove very poorly aligned regions from all the alignments. We will use the `Trim_Alignments_Trimal.py` module to accomplish this. We will use the relabeled fasta alignments found in `/bin/09-Relabel/Output/Relabeled-by-species/Relabeled-fasta-files/`.


We will trim alignments using the following command:

```
python Trim_Alignments_Trimal.py -i /bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/10-Trim/ -f fasta -a gt --gt 0.1
```

This will run the gap-threshold algorithm of trimal, using a threshold value of 0.1. 

The trimmed alignments will be written to the `/bin/10-Trim/Trimmed-fasta-Files/` directory.

#### **Concatenate alignments to create the final supermatrix**

The moment is upon us! It is time to create the final supermatrix. We will use the `Concatenation.py` module to perform this task. We will use the trimmed fasta files for this step, which are in the `/bin/10-Trim/Trimmed-fasta-Files/` directory. 

We will specify the input format is fasta, the output format should be phylip (compatible with RAxML), and that missing characters should be coded with a `-` symbol.

We will perform this task using the following command:

```
python Concatenation.py -i /bin/10-Trim/Trimmed-fasta-Files/ -o /bin/11-Concatenate/ --informat fasta --outformat phylip -s dash
```

And here is the output on screen when we run it:

```
Found 60 fasta files to concatenate.

Found 1,399 unique taxa across alignment files.

Gathering sequences for all taxa (this could take some time)...
	Done.

Writing partitions file.

Writing concatenated alignment.
	Total alignment length = 52,827 bp.
	Total number of sequences included = 12,978.

Finished. Total elapsed time: 0:00:00.620653 (H:M:S)
```

A concatenated alignment in phylip format (`Concatenated_Alignment.phy`), a data partitions file (`Data_Partitions.txt`), and final count of the number of taxa and their number of corresponding loci (`Taxa_Loci_Count.log`) will be written to `/bin/11-Concatenate/`. 

And that's it! We've just created a large supermatrix for iguanian lizards. 

You can use the concatenated alignment to create a phylogeny using RAxML, like we did in our study. Our bootstrapped tree is available in the `OSF/12-RAxML/` folder, in case you'd like to download it.


