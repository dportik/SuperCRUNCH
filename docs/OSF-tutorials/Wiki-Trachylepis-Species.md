## Tutorial for Analysis: Trachylepis-Species

This wiki is intended to allow you to replicate the Trachylepis-Species analysis reported in the Supplemental Materials. All of the input and output files of the completed analysis are provided in the OSF data folders.  

The wiki is presented as a tutorial. However, you can always scroll through the sections of this wiki to see how each module was called, and which files were used and produced at each step. One of the main goals is to provide you with all the input files required for running the different modules of SuperCRUNCH, and to demonstrate the use of various options along the way. This page is therefore intended to function as a tutorial and a reference resource.

**Please note that this analysis requires ~10MB of space and will take ~10 minutes to complete.**

---------------


## **Preparing for the Tutorial**

Before you can run this tutorial locally, there are few directories that must be created, and three files that will need to be downloaded. In total, these files will be ~1MB in size. You will also need ~9MB of additional free space to allow for the analysis files to be created across the tutorial steps.

**The easiest way to prepare is to download the `Directory-Structure-Data.zip` file located in the `OSF/00-Directory-Structure/` folder in the OSF files window. It will contain the complete set of directories and all the data files. Please read the instructions below for more details.** 

### 1. Create Local Directories
To begin, you will need to have the following directories on your local computer:

```
bin/
│
├── 00-Starting-seqs/
├── 01-Get-taxa/
├── 02-Starting-materials/
├── 03-Parse-loci/
├── 04-Cluster-blast/
├── 05-Filter-seqs/
├── 06-Make-acc-table/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Convert/
├── 11-Concatenate/
```

There are two options for getting these folders made. You can download the `Directory-Structure.zip` file or the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder and unzip it. The `Directory-Structure-Data.zip` version contains all three of the required data files already in their correct directories. The `Directory-Structure.zip` version simply contains all the empty directories. Alternatively, if you choose to not download these files you can make the directories by hand. 

Please note that the `bin/` directory is just a placeholder name, to simplify the "path" to each directory. All of the subsequent steps will reference the "paths" to these folders. These paths will be different on your local computer, so make sure to change the path as appropriate. 

In this example, the path to the `02-Starting-materials/` folder will be called `/bin/02-Starting-materials/`. 

### 2. Download Data

You will need to obtain three files from the OSF storage to run the analysis. The easiest way to do this is download the folders with the data already present, by getting the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder. This will have all three files listed below.

Alternatively, if you downloaded the empty directories or made them by hand, you will need to download and move the three files listed below from their OSF directories and move them to their correct local directory.

**Starting sequences:**
The GenBank sequences downloaded in fasta format are available here in the `OSF/00-Fasta-seqs/` folder. The sequences are in the file called `Trachylepis.fasta`. A list of their accession numbers is in the `Accession_List.txt` file. The accession list is not necessary for the analysis.

Download and copy the `Trachylepis.fasta` file to your local `/bin/00-Starting-seqs/` and `/bin/02-Starting-materials/` folder (two locations!). We will use this file in both directories for the analysis.

**Locus search terms and taxon names list:** Download the `Locus-search-terms.txt` and `Taxon-List.txt` files from `OSF/02-Starting-materials/` and move them both to your `/bin/02-Starting-materials/` folder. 

### 3. Check File Locations are Correct

If you downloaded the above files and moved them to the correct directories, you should find them in the following locations. The necessary files are marked with asterisks at the beginning (***) of their names:

```
bin/
│
├── 00-Starting-seqs/
│   └── *** Trachylepis.fasta
├── 01-Get-taxa/
├── 02-Starting-materials/
│   ├── *** Trachylepis.fasta
│   ├── *** Locus-Search-Terms-UCE-5k-set.txt
│   └── *** Taxon-List.txt
├── 03-Parse-loci/
├── 04-Cluster-blast/
├── 05-Filter-seqs/
├── 06-Make-acc-table/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Convert/
├── 11-Concatenate/
```

If everything looks good, you are ready to proceed to the tutorial!

---------------


## **Running the Tutorial**

The goal of this analysis is to create a species-level supermatrix from a published multi-locus phylogeography dataset. This analysis will construct a species-level supermatrix that includes representative sequences for each locus for each species. The other parallel analysis (Trachylepis-Vouchered) can be run to reconstruct the original phylogeographic dataset.

This analysis will take roughly 10 minutes to run (including user-time), and is a good example of how to create a species-level supermatrix quickly. The final supermatrix should contain 7 species, 4 loci, and 26 sequences.

The steps in this tutorial follow along with the steps listed in Supplemental Table S6 of Portik & Wiens (2020), and include:

1. Pre-Analysis - Get Taxon Names
1. Parse Loci
1. Similarity Filtering 
1. Sequence Selection
1. Sequence Alignment
1. Post-Alignment

As a reminder, the paths to folders and files will likely be different on your local computer. You will need to replace the `/bin/` section of the paths used in the commands with the local path to the directories you created in your environment. 

For example, `/bin/02-Starting-materials/` may actually be something like ` /Users/portik/Tutorial-5/02-Starting-materials/`. 

---------------

### **1. Pre-Analysis - Get Taxon Names**

For this analysis, we are targeting < 10 taxa. Since we are reconstructing a published dataset, we could get the taxon names from the paper. However, there is an easier way - taxon names can be extracted directly from fasta files using the `Fasta_Get_Taxa.py` module. 

The fasta file should be located in the `/bin/00-Starting-seqs/` folder. We will use the following command to run the module:

```
python Fasta_Get_Taxa.py -i /bin/00-Starting-seqs/ -o /bin/01-Get-taxa/ 
```

Here is the output on screen:

```
Examining Trachylepis.fasta:

	Read 442 total records...


Found a total of 9 unique 'species' names.
Found a total of 3 unique 'subspecies' names.
```

This will produce two files: `Species_Names.txt` and `Subspecies_Names.txt`. 

I used these two files to create the final `Taxon-List.txt` file in the `/bin/02-Starting-materials/` directory. 

---------------

### **2. Parse Loci**

SuperCRUNCH will first identify loci based on the labels present in the sequence records. All sequences that match a specific gene abbreviation and/or description will be written to a fasta file created for that locus. These searches are conducted using the `Parse_Loci.py` module.

Running this step requires three files: 
- a list of genes/loci to search for
- a list of taxon names
- the fasta sequences to search. 

These files should now all be in the `/bin/02-Starting-materials/` directory, and include `Locus-search-terms.txt`, `Taxon-List.txt`, and `Trachylepis.fasta`. 

We will run `Parse_Loci.py` with the following command:

```
python Parse_Loci.py -i /bin/02-Starting-materials/Trachylepis.fasta -l /bin/02-Starting-materials/Locus-search-terms.txt -t /bin/02-Starting-materials/Taxon-List.txt -o /bin/03-Parse-loci/ 
```

For this analysis, we are not interested in finding subspecies. If the taxon list file contains subspecies names, we would use the `--no_subspecies` flag. However, the taxon list file used here only contains species names, and so whether or not we use the `--no_subspecies` flag, the outcome will be the same.

All outputs are written to the `/bin/03-Parse-loci/` folder.

Inside there are two new folders:
- `Summary-File/` contains a text file showing the number of sequences found per locus.
- `Parsed-Fasta-Files/` contains all the locus-specific fasta files for which at least two sequences were found.

The SQL database generated during this step will also be in the `/bin/03-Parse-loci/` folder.

The locus-specific fasta files contained in the `/bin/03-Parse-loci/Parsed-Fasta-Files/` folder will be used for the next step.

---------------

### **3. Similarity Filtering**

SuperCRUNCH performs similarity searches with BLASTn. It offers two options, which differ in whether the reference sequences are selected automatically or provided by the user. 

The automatic reference option is appropriate for loci containing "simple" records. We define “simple” record sets as those generally containing a single gene region with limited length variation, which results from use of the same primers (Sanger-sequencing) or probes (sequence capture) to generate sequences. 

The user-supplied reference option is appropriate for loci containing "complex" records. We define “complex” records as those containing the target region plus non-target sequence (e.g., other regions or genes). Common examples include long mtDNA fragments and whole mitogenomes, and genes that have been sequenced for different fragments that have little or no overlap. 

Sequences that have significant hits to the references will be trimmed if necessary. This means that only the target regions of the input sequences are included in the output fasta file. 

#### **Run Cluster_Blast_Extract**

Similarity searches using automatically selected references can be performed with the `Cluster_Blast_Extract.py` module. 

For this analysis, nearly all of the sequences we downloaded came from a single study. This means that those sequences were generated using the same methods, and as a result they should be fairly uniform in length. For this reason, they will fall into the "simple" records category. That is why we will use this particular filtering strategy for all 4 loci (including the mtDNA locus ND2).

All four locus-specific fasta files will be located in the `/bin/03-Parse-loci/Parsed-Fasta-Files/` folder.

We will run the module with the following command:

```
python Cluster_Blast_Extract.py -i /bin/03-Parse-loci/Parsed-Fasta-Files/ -o /bin/04-Cluster-blast/ -b dc-megablast -m span --threads 4
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Four new directories have been written to the `/bin/04-Cluster-blast/` folder, which contain all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/04-Cluster-blast/Filtered-Fasta-Files/` directory. We will use these fasta files for the sequence selection step. 

---------------

### **4. Sequence Selection**


Here we select sequences to create a species-level supermatrix, and make a table of GenBank accession numbers.

#### **Select Sequences**

For supermatrices, a single sequence is typically used to represent each species for each gene. If multiple sequences for a given gene exist for a given species (e.g., because multiple individuals were sampled), then an objective strategy must be used for sequence selection. We can use the `Filter_Seqs_and_Species.py` module to accomplish this. 

We will select sequences for all loci based only on length (longest). We will also set 200 bp as the minimum length required to keep a sequence.

To run `Filter_Seqs_and_Species.py`, we need the fasta files that passed the similarity filtering (which are located in `/bin/04-Cluster-blast/Filtered-Fasta-Files/`), and a taxon names list file (`Taxon-List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/04-Cluster-blast/Filtered-Fasta-Files/ -s oneseq -f length -m 200 -t /bin/02-Starting-materials/Taxon-List.txt -o /bin/05-Filter-seqs/ --quiet
```

Notice that for this analysis, we exclude subspecies names by using the `--no_subspecies` flag.

The outputs are written to three folders contained in `/bin/05-Filter-seqs/Results/`. 

The main outputs are the fasta files with the filtered sequences, which are located in `/bin/05-Filter-seqs/Results/Filtered-Fasta-Files/`. We will use these for the next steps.


#### **Make a table of GenBank accession numbers**

Organizing GenBank accession numbers can be a serious undertaking, but this can be automated for you by using the `Make_Acc_Table.py` module. 

This step can be performed any time after sequence selection, and here we perform it after the similarity filtering step.

We can generate this handy table using the following command:

```
python Make_Acc_Table.py -i /bin/05-Filter-seqs/Results/Filtered-Fasta-Files/ -o /bin/06-Acc-table/
```

By default, subspecies names are not included.

The outputs in `/bin/06-Make-acc-table/` include the table itself (`GenBank_Accession_Table.txt`) and a summary of the sample included and the number of sequences/loci for each sample (`Taxon_Summary.txt`).

---------------

### **5. Sequence Alignment**

In preparation for alignment, we need to complete one pre-alignment step. We will adjust the direction of all sequences within each fasta file. This will ensure they are all oriented in the same direction, and prevent catastrophes in sequence alignment.


#### **Pre-Alignment: Adjust sequence directions**

We will ensure the sequences in each fasta file are correctly oriented by using the `Adjust_Direction.py` module. All of the target fasta files should be in the same output directory from sequence selection (`//bin/05-Filter-seqs/Results/Filtered-Fasta-Files/`), so we will just run this step using that location. 

Here is the command we will use:

```
python Adjust_Direction.py -i /bin/05-Filter-seqs/Results/Filtered-Fasta-Files/ -o /bin/07-Adjust/ --threads 8 
```

You should adjust the number of threads accordingly.

There should now be two output directories and a summary file in the `/bin/07-Adjust/` folder. The adjusted fasta files are the main output we are interested in, and they are located in the `/bin/07-Adjust/Adjusted-Fasta-Files/` folder. We will use these for alignment.

#### **Multiple Sequence Alignment**

We will use MAFFT (--auto) to align all loci. This is by far the fastest alignment option, and should take <20 seconds total (with 8 threads!).

We will perform sequence alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/07-Adjust/Adjusted-Fasta-Files/ -o /bin/08-Align/ -a mafft --threads 8 
```

You should adjust the number of threads accordingly.

The final alignments will be written to `/bin/08-Align/Alignments-MAFFT/`.

---------------

### **6. Post-Alignment**

We will now perform three post-alignment tasks:

- We will relabel the sequence records with relevant taxon names, making them compatible with downstream programs.
- We will convert the fasta alignments into nexus and phylip format.
- We will concatenate the alignments to make a final supermatrix.

#### **Relabel sequences in alignments**

Currently, the alignments contain sequences that are labeled using a modified version of the GenBank record line. This cannot be used for most downstream programs. We will relabel these sequences to include the taxon label.

All of the target fasta files are in the `/bin/08-Align/Alignments-MAFFT/` folder.

We will relabel sequences with the `Fasta_Relabel_Seqs.py` module using the following command:

```
python Fasta_Relabel_Seqs.py -i /bin/08-Align/Alignments-MAFFT/ -o /bin/09-Relabel/ -r species
```

After this step, sequences will be labeled by a taxon name (such as `Trachylepis_punctulata`). We did not allow for subspecies names here. 

The fasta alignments with relabeled sequences will be written to `/bin/09-Relabel/Relabeled-by-species/`.


#### **Convert alignments into nexus and phylip format**

We can quickly convert any relabeled fasta alignment to nexus and phylip format using the `Fasta_Convert.py` module. We will perform this format conversion for all of the relabeled alignments found in `/bin/09-Relabel/Relabeled-by-species/`.

We will convert formats using the following command:

```
python Fasta_Convert.py -i /bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/10-Convert/
```

The nexus-formatted files will be written to `/bin/10-Convert/Nexus-Files/` and the phylip-formatted files will be written to `/bin/10-Convert/Phylip-Files/`.



#### **Concatenate alignments to create the final supermatrix**

The moment is upon us! It is time to create the final species-level supermatrix. We will use the `Concatenation.py` module to perform this task. We will use the relabeled fasta files for this step, which are in the `/bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/` directory. 

We will specify the input format is fasta, the output format should be phylip (compatible with RAxML), and that missing characters should be coded with a `-` symbol.

We will perform this task using the following command:

```
python Concatenation.py -i /bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/11-Concatenate/ --informat fasta --outformat phylip -s dash
```

And here is the output on screen when we run it:

```
Found 4 fasta files to concatenate.

Found 7 unique taxa across alignment files.

Gathering sequences for all taxa (this could take some time)...
	Done.

Writing partitions file.

Writing concatenated alignment.
	Total alignment length = 3,697 bp.
	Total number of sequences included = 26.

Finished. Total elapsed time: 0:00:00.006848 (H:M:S)
```

A concatenated alignment in phylip format (`Concatenated_Alignment.phy`), a data partitions file (`Data_Partitions.txt`), and final count of the number of taxa and their number of corresponding loci (`Taxa_Loci_Count.log`) will be written to `/bin/10-Concatenate/`. 

And that's it! 

This is a very small supermatrix (7 species, 4 loci). It is primarily intended to illustrate a major difference in sequence selection for the same dataset. Compare these results to those obtained in the Trachylepis-Vouchered tutorial. The only difference between the analyses was selecting one representative sequence per species per locus, versus selecting all vouchered sequences per species per locus.

You can use the concatenated alignment to create a phylogeny using RAxML, like we did in our study. Our bootstrapped tree is available in the `OSF/11-RAxML/` folder, in case you'd like to download it.


