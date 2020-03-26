## Tutorial for Analysis: Kaloula-Species

This wiki is intended to allow you to replicate the Kaloula-Species analysis reported in the Supplemental Materials. All of the input and output files of the completed analysis are provided in the OSF data folders.  

The wiki is presented as a tutorial. However, you can always scroll through the sections of this wiki to see how each module was called, and which files were used and produced at each step. One of the main goals is to provide you with all the input files required for running the different modules of SuperCRUNCH, and to demonstrate the use of various options along the way. This page is therefore intended to function as a tutorial and a reference resource.

**Please note that this analysis requires ~300MB of space and will take at least 20 minutes to complete.**

---------------


## **Preparing for the Tutorial**

Before you can run this tutorial locally, there are few directories that must be created, and three files that will need to be downloaded. In total, these files will be ~35MB in size. You will also need ~250MB of additional free space to allow for the analysis files to be created across the tutorial steps.

**The easiest way to prepare is to download the `Directory-Structure-Data.zip` file located in the `OSF/00-Directory-Structure/` folder in the OSF files window. It will contain the complete set of directories and all the data files. Please read the instructions below for more details.** 

### 1. Create Local Directories
To begin, you will need to have the following directories on your local computer:

```
bin/
│
├── 00-Fasta-seqs/
├── 01-Get-taxa/
├── 02-Starting-materials/
├── 03-Parse-loci/
├── 04-Filter-seqs/
├── 05-Make-acc-table/
├── 06-Adjust/
├── 07-Align/
├── 08-Relabel/
├── 09-Concatenate/
```

There are two options for getting these folders made. You can download the `Directory-Structure.zip` file or the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder and unzip it. The `Directory-Structure-Data.zip` version contains all three of the required data files already in their correct directories. The `Directory-Structure.zip` version simply contains all the empty directories. Alternatively, if you choose to not download these files you can make the directories by hand. 

Please note that the `bin/` directory is just a placeholder name, to simplify the "path" to each directory. All of the subsequent steps will reference the "paths" to these folders. These paths will be different on your local computer, so make sure to change the path as appropriate. 

In this example, the path to the `02-Starting-materials/` folder will be called `/bin/02-Starting-materials/`. 

### 2. Download Data

You will need to obtain three files from the OSF storage to run the analysis. The easiest way to do this is download the folders with the data already present, by getting the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder. This will have all three files listed below.

Alternatively, if you downloaded the empty directories or made them by hand, you will need to download and move the three files listed below from their OSF directories and move them to their correct local directory.

**Starting sequences:**
The GenBank sequences downloaded in fasta format are available here in the `OSF/00-Fasta-seqs/` folder. The sequences are in the file called `Anura_UCE.fasta`. A list of their accession numbers is in the `Accession_List.txt` file. The accession list is not necessary for the analysis.

Download and copy the `Anura_UCE.fasta` file to your local `/bin/00-Fasta-seqs/` and `/bin/02-Starting-materials/` folder (two locations!). We will use this file in both directories for the analysis.

**Locus search terms and taxon names list:** Download the `Locus-Search-Terms-UCE-5k-set.txt` and `Taxon-List.txt` files from `OSF/02-Starting-materials/` and move them both to your `/bin/02-Starting-materials/` folder. 

### 3. Check File Locations are Correct

If you downloaded the above files and moved them to the correct directories, you should find them in the following locations. The necessary files are marked with asterisks at the beginning (***) of their names:

```
bin/
│
├── 00-Fasta-seqs/
│   └── *** Anura_UCE.fasta
├── 01-Get-taxa/
├── 02-Starting-materials/
│   ├── *** Anura_UCE.fasta
│   ├── *** Locus-Search-Terms-UCE-5k-set.txt
│   └── *** Taxon-List.txt
├── 03-Parse-loci/
├── 04-Filter-seqs/
├── 05-Make-acc-table/
├── 06-Adjust/
├── 07-Align/
├── 08-Relabel/
├── 09-Concatenate/
```

If everything looks good, you are ready to proceed to the tutorial!

---------------


## **Running the Tutorial**

The goal of this analysis is to reconstruct a phylogenomic supermatrix from published UCE data. The original dataset contained a mix of species and population-level sampling. This analysis will construct a species-level UCE supermatrix that only includes a single representative for each species/subspecies. The other parallel analysis (Kaloula-Vouchered) is intended to include all of the samples.

This analysis will take roughly 20 minutes to run (not including user-time), and is a good example of how to create a phylogenomic supermatrix quickly. The final supermatrix should contain 14 taxa, 1,784 UCE loci, and 22,717 sequences.

The steps in this tutorial follow along with the steps listed in Supplemental Table S4 of Portik & Wiens (2020), and include:

1. Pre-Analysis - Get Taxon Names
1. Parse Loci
1. Sequence Selection
1. Sequence Alignment
1. Post-Alignment

As a reminder, the paths to folders and files will likely be different on your local computer. You will need to replace the `/bin/` section of the paths used in the commands with the local path to the directories you created in your environment. 

For example, `/bin/02-Starting-materials/` may actually be something like ` /Users/portik/Tutorial-Iguania/02-Starting-materials/`. 

---------------

### **1. Pre-Analysis - Get Taxon Names**

For this analysis, we are targeting < 20 taxa. Since we are reconstructing a published dataset, we could get the taxon names from the paper. However, there is an easier way - taxon names can be extracted directly from fasta files using the `Fasta_Get_Taxa.py` module. 

The fasta file should be located in the `/bin/00-Fasta-seqs/` folder. We will use the following command to run the module:

```
python Fasta_Get_Taxa.py -i /bin/00-Fasta-seqs/ -o /bin/01-Get-taxa/ 
```


Here is the output on screen:

```
Examining Anura_UCE.fasta:
	Processed 10,000 records...
	Processed 20,000 records...
	Processed 30,000 records...

	Read 38,568 total records...


Found a total of 12 unique 'species' names.
Found a total of 4 unique 'subspecies' names.
```

This will produce two files: `Species_Names.txt` and `Subspecies_Names.txt`. 

I used these two files to create the final `Taxon-list.txt` file in the `/bin/02-Starting-materials/` directory. 

---------------

### **2. Parse Loci**

SuperCRUNCH will first identify loci based on the labels present in the sequence records. All sequences that match a specific gene abbreviation and/or description will be written to a fasta file created for that locus. These searches are conducted using the `Parse_Loci.py` module.

Running this step requires three files: 
- a list of genes/loci to search for
- a list of taxon names
- the fasta sequences to search. 

These files should now all be in the `/bin/02-Starting-materials/` directory, and include `Locus-Search-Terms-UCE-5k-set.txt`, `Taxon-List.txt`, and `Anura_UCE.fasta`. 

We will run `Parse_Loci.py` with the following command:

```
python Parse_Loci.py -i /bin/02-Starting-materials/Anura_UCE.fasta -l /bin/02-Starting-materials/Locus-Search-Terms-UCE-5k-set.txt -t /bin/02-Starting-materials/Taxon-List.txt -o /bin/03-Parse-loci/ 
```
Notice that for this analysis, we want to search for subspecies, and consequently we did not use the `--no_subspecies` flag.

All outputs are written to the `/bin/03-Parse-loci/` folder.

Inside there are two new folders:
- `Summary-File/` contains a text file showing the number of sequences found per locus.
- `Parsed-Fasta-Files/` contains all the locus-specific fasta files for which at least two sequences were found.

The SQL database generated during this step will also be in the `/bin/03-Parse-loci/` folder.

The locus-specific fasta files contained in the `/bin/03-Parse-loci/Parsed-Fasta-Files/` folder will be used for the next step.

---------------

**PLEASE NOTE - In general the locus-parsing step should always be followed by a similarity searching filter (via BLAST).**

**This helps ensure the sequences are actually orthologous, and also trims off non-target regions.**

**However, because these UCE sequences already passed through a similarity filter using PHYLUCE, we will skip this step in SuperCRUNCH.**

---------------

### **3. Sequence Selection**


Here we select a single sequence per species per locus, and make a table of GenBank accession numbers.

#### **Select Representative Sequences**

For supermatrices, a single sequence is typically used to represent each species for each gene. If multiple sequences exist for a gene for a given species/subspecies (e.g., because multiple individuals were sampled), then an objective strategy must be used for sequence selection. We can use the `Filter_Seqs_and_Species.py` module to accomplish this. 

We will select sequences for all loci based only on length (longest). We will also set 150 bp as the minimum length required to keep a sequence.

To run `Filter_Seqs_and_Species.py`, we need the locus-specific fasta files (which are located in `/bin/03-Parse-loci/Parsed-Fasta-Files/`), and a taxon names list file (`Taxon-List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/03-Parse-loci/Parsed-Fasta-Files/ -o /bin/04-Filter-seqs/ -s oneseq -f length -m 150 -t /bin/02-Starting-materials/Taxon-List.txt --quiet
```

Notice that for this analysis, we want to search for subspecies, and consequently we did not use the `--no_subspecies` flag.

The outputs are written to three folders contained in `/bin/04-Filter-seqs/Results/`. 

The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/04-Filter-seqs/Results/Filtered-Fasta-Files/`. We will use these for the next steps.


#### **Make a table of GenBank accession numbers**

Organizing GenBank accession numbers can be a serious undertaking, but this can be automated for you by using the `Make_Acc_Table.py` module. 

This step can be performed any time after sequence selection, and here we perform it right after the sequence selection step.

We can generate this handy table using the following command:

```
python Make_Acc_Table.py -i /bin/04-Filter-seqs/Results/Filtered-Fasta-Files/ -o /bin/05-Make-acc-table/ -s /bin/02-Starting-materials/Taxon-List.txt 
```

We included subspecies in this table using the `-s` flag, which requires the taxon list so the subspecies names can be detected. We can confirm subspecies are included by looking at the output on screen:

```
Using subspecies labels found in: Taxon-List.txt

Getting subspecies from Taxon-List.txt.
	Found 4 subspecies names to include.

Found 1,784 files to use for accession table.

Found 14 unique taxon labels across fasta files.

Gathering accession numbers for all taxa.

Writing accession table.

Writing taxon summary file.
```

The outputs in `/bin/05-Make-acc-table/` include the table itself (`GenBank_Accession_Table.txt`) and a summary of the taxa included and the number of sequences/loci for each taxon (`Taxon_Summary.txt`).

Given that this table includes ~1,800 loci, it would have been extremely difficult to organize these accession numbers any other way. 

---------------

### **4. Sequence Alignment**

In preparation for alignment, we need to complete one pre-alignment step. We will adjust the direction of all sequences within each fasta file. This will ensure they are all oriented in the same direction, and prevent catastrophes in sequence alignment.


#### **Pre-Alignment: Adjust sequence directions**

We will ensure the sequences in each fasta file are correctly oriented by using the `Adjust_Direction.py` module. All of the target fasta files should be in the same output directory from sequence selection (`/bin/04-Filter-seqs/Results/Filtered-Fasta-Files/`), so we will just run this step using that location. 

Here is the command we will use:

```
python Adjust_Direction.py -i /bin/04-Filter-seqs/Results/Filtered-Fasta-Files/ -o /bin/06-Adjust/ --threads 8
```

You should adjust the number of threads accordingly.

There should now be two output directories and a summary file in the `/bin/06-Adjust/` folder. The adjusted fasta files are the main output we are interested in, and they are located in the `/bin/06-Adjust/Adjusted-Fasta-Files/` folder. We will use these for alignment.

#### **Multiple Sequence Alignment**

We will use MAFFT (--auto) to align all loci. This is by far the fastest alignment option, and should take ~15 minutes total (with 8 threads!).

We will perform sequence alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/06-Adjust/Adjusted-Fasta-Files/ -o /bin/07-Align/ -a mafft --threads 8
```

You should adjust the number of threads accordingly.

The final alignments will be written to `/bin/07-Align/Alignments-MAFFT/`.

---------------

### **5. Post-Alignment**

We will now perform two post-alignment tasks:

- We will relabel the sequence records with relevant taxon names, making them compatible with downstream programs.
- We will concatenate the alignments to make a final supermatrix.

#### **Relabel sequences in alignments**

Currently, the alignments contain sequences that are labeled using a modified version of the GenBank record line. This cannot be used for most downstream programs. We will relabel these sequences to include only the taxon label, and will include subspecies names.

All of the target fasta files are in the `/bin/07-Align/Alignments-MAFFT/` folder.

We will relabel sequences with the `Fasta_Relabel_Seqs.py` module using the following command:

```
python Fasta_Relabel_Seqs.py -i /bin/07-Align/Alignments-MAFFT/ -o /bin/08-Relabel/ -r species -s /bin/02-Starting-materials/Taxon-List.txt 
```

Notice that to include subspecies names, we must again provide the taxon list file (the one used for previous steps). 

The fasta alignments with relabeled sequences will be written to `/bin/08-Relabel/Relabeled-by-species/`.


#### **Concatenate alignments to create the final supermatrix**

The moment is upon us! It is time to create the final supermatrix. We will use the `Concatenation.py` module to perform this task. We will use the relabeled fasta files for this step, which are in the `/bin/08-Relabel/Relabeled-by-species/` directory. 

We will specify the input format is fasta, the output format should be phylip (compatible with RAxML), and that missing characters should be coded with a `-` symbol.

We will perform this task using the following command:

```
python Concatenation.py -i /bin/08-Relabel/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/09-Concatenate/ --informat fasta --outformat phylip -s dash
```

And here is the output on screen when we run it:

```
Found 1,784 fasta files to concatenate.

Found 14 unique taxa across alignment files.

Gathering sequences for all taxa (this could take some time)...
	Done.

Writing partitions file.

Writing concatenated alignment.
	Total alignment length = 1,352,535 bp.
	Total number of sequences included = 22,717.

Finished. Total elapsed time: 0:00:00.809911 (H:M:S)
```

A concatenated alignment in phylip format (`Concatenated_Alignment.phy`), a data partitions file (`Data_Partitions.txt`), and final count of the number of taxa and their number of corresponding loci (`Taxa_Loci_Count.log`) will be written to `/bin/09-Concatenate/`. 

And that's it! We've just created reconstructed a ~1,800 locus UCE supermatrix. 

You can use the concatenated alignment to create a phylogeny using RAxML, like we did in our study. Our bootstrapped tree is available in the `OSF/10-RAxML/` folder, in case you'd like to download it.


