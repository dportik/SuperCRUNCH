## Tutorial for Analysis: Hyperoliid-Outgroup 

This wiki is intended to allow you to replicate the Hyperoliid-Outgroup analysis reported in the Supplemental Materials. All of the input and output files of the completed analysis are provided in the OSF data folders.  

The wiki is presented as a tutorial. However, you can always scroll through the sections of this wiki to see how each module was called, and which files were used and produced at each step. One of the main goals is to provide you with all the input files required for running the different modules of SuperCRUNCH, and to demonstrate the use of various options along the way. This page is therefore intended to function as a tutorial and a reference resource.

**Please note that this analysis requires ~100MB of space and will take at least ~10 minutes to complete.**

---------------


## **Preparing for the Tutorial**

Before you can run this tutorial locally, there are several directories that must be created, and a few files that will need to be downloaded. In total, these files will be ~10MB in size. You will also need ~90MB of additional free space to allow for the analysis files to be created across the tutorial steps.

**The easiest way to prepare is to download the `Directory-Structure-Data.zip` file located in the `OSF/00-Directory-Structure/` folder in the OSF files window. It will contain the complete set of directories and all the data files (except one, see below!). Please read the instructions below for more details.** 

### 1. Create Local Directories
To begin, you will need to have the following directories on your local computer:

```
bin/
│
├── 01-Starting-seqs/
│	├── GenBank/
│	└── Local/
├── 02-Get-taxa/
│	├── Arthroleptidae/
│	└── Hyperoliidae/
├── 03-Starting-materials/
├── 04-Parse-loci/
│	├── All/
│	└── Hyperoliids/
├── 05-Reference-blast/
│	├── Input/
│	└── Output/
├── 06-Filter-seqs/
│	├── Input/
│	└── Output/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Concatenate/
```
There are two options for getting these folders made. You can download the `Directory-Structure.zip` file or the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder and unzip it. The `Directory-Structure-Data.zip` version contains all but one of the required data files already in their correct directories (see below). The `Directory-Structure.zip` version simply contains all the empty directories. Alternatively, if you choose to not download these files you can make the directories by hand. 

Please note that the `bin/` directory is just a placeholder name, to simplify the "path" to each directory. All of the subsequent steps will reference the "paths" to these folders. These paths will be different on your local computer, so make sure to change the path as appropriate. 

In this example, the path to the `01-Starting-materials/` folder will be called `/bin/01-Starting-materials/`. 

### 2. Download Data

You will need to obtain six files from the OSF storage to run the analysis. The easiest way to do this is download the folders with the data already present, by getting the `Directory-Structure-Data.zip` file from the `OSF/00-Directory-Structure/` folder. This will have all of the files listed below, except for the starting sequences. You will need to download the starting sequences before running the tutorial!

Alternatively, if you downloaded the empty directories or made them by hand, you will need to download and move the files listed below from their OSF directories and move them to their correct local directory.

**Starting sequences:**
There are three different sets of sequences we will use: `Arthroleptidae.fasta`, `Hyperoliidae.fasta`, and `Arthroleptidae-Hyperoliidae.fasta`. 
 - Download `Arthroleptidae.fasta` from `OSF/01-Starting-Seqs/GenBank/` and move it to your local `/bin/01-Starting-Seqs/GenBank/` folder.
 - Download `Hyperoliidae.fasta` from `OSF/01-Starting-Seqs/Local/` and move it to your local `/bin/01-Starting-Seqs/Local/` folder.
 - Download `Arthroleptidae-Hyperoliidae.fasta` from `OSF/03-Starting-materials/` and move it to your local `/bin/01-Starting-materials/` folder.


**Locus search terms and taxon names list:** Download the `Locus-search-terms.txt` and `Taxon-List.txt` files from `OSF/03-Starting-materials/` and move them both to your `/bin/03-Starting-materials/` folder. 

**Files for Reference Blast:** Download the `multisearch-key.txt` file from `OSF/05-Reference-blast/` and move it to your `/bin/05-Reference-blast/` folder.

### 3. Check File Locations are Correct

If you downloaded the above files and moved them to the correct directories, you should find them in the following locations. The downloaded files are marked with asterisks at the beginning (***) of their names:

```
bin/
│
├── 01-Starting-seqs/
│	├── GenBank/
│   |   └── *** Arthroleptidae.fasta
│	└── Local/
│       └── *** Hyperoliidae.fasta
├── 02-Get-taxa/
│	├── Arthroleptidae/
│	└── Hyperoliidae/
├── 03-Starting-materials/
│   ├── *** Arthroleptidae-Hyperoliidae.fasta
│   ├── *** Locus-search-terms.txt
│   └── *** Taxon-List.txt
├── 04-Parse-loci/
│	├── All/
│	└── Hyperoliids/
├── 05-Reference-blast/
│	├── *** multisearch-key.txt
│	├── Input/
│	└── Output/
├── 06-Filter-seqs/
│	├── Input/
│	└── Output/
├── 07-Adjust/
├── 08-Align/
├── 09-Relabel/
├── 10-Concatenate/
```

If everything looks good, you are ready to proceed to the tutorial.

---------------


## **Running the Tutorial**

Here we use SuperCRUNCH to perform a common but sometimes exceedingly difficult task in phylogenetics: adding published outgroup sequences to an unpublished sequencing project. The local dataset consists of six loci sequenced for ~128 species, but many species are represented by multiple vouchered samples. For this analysis, we want to add all available GenBank sequence data for the family Arthroleptidae, which is the sister family of Hyperoliidae (the ingroup). For the local sequences, we want to treat all samples as distinct (equivalent to a vouchered analysis), whereas for the outgroups we simply want to include all possible data for a given species (equivalent to a species-level analysis, in which sequences for a species most likely come from different samples). 

Within the local sequences, I intentionally labeled the records such that the species names was followed immediately by a museum/field identifier, which allows us to take advantage of the flexible “subspecies” option in SuperCRUNCH. For information on how to label your local sequences please see [**here**](https://github.com/dportik/SuperCRUNCH/wiki/2:-Starting-Sequences#ULS). To label my sequences, I wrote a custom Python script to help. You may do the same, or choose to use something like Geneious to help format the sequence labels. For information about how taxon labels can be searched for (including by using voucher codes as part of the name), please see [**here**](https://github.com/dportik/SuperCRUNCH/wiki/3:-Assessing-Taxonomy). The subspecies option allows any three-part name to be used, and the third part of a name can contain either a subspecies label or any alphanumerical identifier. 

One interesting strategy we will implement is to use all of the sequences of the hyperoliids (the unpublished data) to act as references during similarity searches of the combined sequences (arthroleptids and hyperoliids). This will ensure that the arthroleptid sequences match the target regions of the loci for the hyperoliids. 

In the end, we will create a hybrid supermatrix - the ingroup can have multiple representatives per species (labeled with voucher codes), whereas the outgroup can only have a single representative per species. The final hybrid matrix will include 6 loci, 365 samples, and 1,724 sequences. 

This analysis will take roughly ~10 minutes to run, but you will likely spend extra time preparing for certain steps. 

This analysis is a good example of how to merge published and unpublished sequences to complete a large phylogenetics project.

The steps in this tutorial follow along with the steps listed in Supplemental Table S9 of Portik & Wiens (2020), and include:

1. Pre-Analysis - Get Taxon Names
1. Parse Loci
1. Similarity Filtering
1. Sequence Selection
1. Sequence Alignment
1. Post-Alignment

As a reminder, the paths to folders and files will likely be different on your local computer. You will need to replace the `/bin/` section of the paths used in the commands with the local path to the directories you created in your environment. 

For example, `/bin/03-Starting-materials/` may actually be something like `/Users/portik/Tutorial-9/03-Starting-materials/`. 

---------------

### **1. Pre-Analysis - Get Taxon Names**

For this analysis, we are dealing with two fundamentally different sequence sets.

`Arthroleptidae.fasta` contains GenBank sequences, and we want species-level data from this sequence set. To begin, we will create a list of the species contained in `Arthroleptidae.fasta`.

`Hyperoliidae.fasta` contains unpublished sequences, and we want to include ALL available sequences for each species. Sequences belonging to the same species are distinguishable based on voucher codes, which we will include as part of their name. SuperCRUNCH has a feature that allows the taxon + voucher code to be considered a three-part name. This is analagous to a "subspecies" label, but the voucher code is in the third part of the name, rather than a scientific name. For more information on this topic, please see [**here**](https://github.com/dportik/SuperCRUNCH/wiki/3:-Assessing-Taxonomy). We want to extract a list of these three-part names (e.g., genus species voucher) from the `Hyperoliidae.fasta` fasta file.

Since we are targeting two different files with different desired outcomes, we will run the `Fasta_Get_Taxa.py` module twice. 

To extract species names from `Arthroleptidae.fasta`, we will use the following command:

```
python Fasta_Get_Taxa.py -i /bin/01-Starting-seqs/GenBank/Arthroleptidae.fasta -o /bin/02-Get-taxa/Arthroleptidae/ 
```

To extract the three-part names (that allow for voucher codes) from `Hyperoliidae.fasta`, we will use the following command:

```
python Fasta_Get_Taxa.py -i /bin/01-Starting-seqs/Local/Hyperoliidae.fasta -o /bin/02-Get-taxa/Hyperoliidae/ --numerical
```

The `--numerical` flag will allow for voucher codes in the third part of the three-part names created. Normally, it will exclude voucher codes in an effort to generate "sensible" subspecies labels, which should not have numerical characters in their name!

For both implementations, two files will be written to the relevant output folder: `Species_Names.txt` and `Subspecies_Names.txt`. 

I used all the three-part names generate for the hyperoliids, and all the species names of the arthroleptids, to make the final `Taxon-List.txt` file that is in the `/bin/03-Starting-materials` directory. I encourage you to look at this list to understand what we will be targeting in the next steps. This strategy will allow us to incorporate ALL the hyperoliid samples while targeting species-level sampling for the arthroleptids.

---------------

### **2. Parse Loci**

SuperCRUNCH will first identify loci based on the labels present in the sequence records. All sequences that match a specific gene abbreviation and/or description will be written to a fasta file created for that locus. These searches are conducted using the `Parse_Loci.py` module.

Running this step requires three files: 
- a list of genes/loci to search for
- a list of taxon names
- the fasta sequences to search. 

Our overall strategy is to use the hyperoliid sequences as references for the similarity searches of the full dataset (hyperoliids + arthroleptids). At this step, we will generate these reference sets, and generate the loci from the full dataset.

We will use the `Parse_Loci.py` module twice:
1. We will create locus-specific fasta files using ONLY the hyperoliid (unpublished) sequences in `Hyperoliidae.fasta`. These will serve as the reference sequence sets downstream.
2. We will create locus-specific fasta files using ALL the sequences (arthroleptids + hyperoliids) in `Arthroleptidae-Hyperoliidae.fasta`. These are the sequences that we will screen against the references created in step 1.


We will run `Parse_Loci.py` to create the reference sequence sets from hyperoliids with the following command:

```
python Parse_Loci.py -i /bin/01-Starting-seqs/Local/Hyperoliidae.fasta -l /bin/03-Starting-materials/Locus-search-terms.txt -t /bin/03-Starting-materials/Taxon-List.txt -o /bin/04-Parse-loci/Hyperoliids/ 
```

We will run `Parse_Loci.py` to create the locus-specific fasta files from the combined sequence sets using the following command:

```
python Parse_Loci.py -i /bin/03-Starting-materials/Arthroleptidae-Hyperoliidae.fasta -l /bin/03-Starting-materials/Locus-search-terms.txt -t /bin/03-Starting-materials/Taxon-List.txt -o /bin/04-Parse-loci/All/
```

This will provide us with all the materials we need to perform the similarity filtering step.

---------------

### **3. Similarity Filtering**

SuperCRUNCH performs similarity searches with BLASTn. It offers two options, which differ in whether the reference sequences are selected automatically or provided by the user. 

The automatic reference option is appropriate for loci containing "simple" records. We define “simple” record sets as those generally containing a single gene region with limited length variation, which results from use of the same primers (Sanger-sequencing) or probes (sequence capture) to generate sequences. 

The user-supplied reference option is appropriate for loci containing "complex" records. We define “complex” records as those containing the target region plus non-target sequence (e.g., other regions or genes). Common examples include long mtDNA fragments and whole mitogenomes, and genes that have been sequenced for different fragments that have little or no overlap. 

Sequences that have significant hits to the references will be trimmed if necessary. This means that only the target regions of the input sequences are included in the output fasta file. 


#### **Run Reference_Blast_Extract using Hyperoliid sequences as references**

We will use the hyperoliid sequences as references to perform similarity searches. This will ensure the all of the published sequence data (Arthroleptids) will be trimmed to the target regions present in the hyperoliid sequences. Similarity searches using user-supplied references can be performed with the `Reference_Blast_Extract.py` module. 

The hyperoliid fasta files we will use as reference sets are contained in the `/bin/04-Parse-loci/Hyperoliids/Parsed-Fasta-Files/` folder. **You will need to relabel these files so that we can use them as reference sequence sets.** For each of the 6 output fasta files in this folder, simply add `H-` to the front of the name. For example, relabel `16S.fasta` to `H-16S.fasta`. Do this for each fasta file. Then copy or move these files to the `/bin/05-Reference-blast/Input/` folder. 

The fasta files we will screen against the references can be found in the `/bin/04-Parse-loci/All/Parsed-Fasta-Files/` folder. Copy or move all 6 of these fasta files to the `/bin/05-Reference-blast/Input/` folder.

This should result in 12 files present in the `/bin/05-Reference-blast/Input/` folder. Just to confirm, this is what should now be in the directories:

```
bin/
│
├── 05-Reference-blast/
│	├── *** multisearch-key.txt
│	├── Input/
│   |   ├── *** H-16S.fasta
│   |   ├── *** H-FICD.fasta
│   |   ├── *** H-KIAA2013.fasta
│   |   ├── *** H-POMC.fasta
│   |   ├── *** H-RAG1.fasta
│   |   ├── *** H-TYR.fasta
│   |   ├── *** 16S.fasta
│   |   ├── *** FICD.fasta
│   |   ├── *** KIAA2013.fasta
│   |   ├── *** POMC.fasta
│   |   ├── *** RAG1.fasta
│   |   └── *** TYR.fasta
│	└── Output/
```

We will take advantage of the multisearch option of `Reference_Blast_Extract.py` to automate runs. This will pair each fasta file with a reference file and run them sequentially. To use this option we will need the `multisearch-key.txt` file, which should be located in the `/bin/05-Reference-blast/` folder.

We are now ready to run the `Reference_Blast_Extract.py` module. 

We will run the module with the following command:

```
python Reference_Blast_Extract.py -i /bin/05-Reference-blast/Input/ -o /bin/05-Reference-blast/Output/ --multisearch /bin/05-Reference-blast/multisearch-key.txt -b dc-megablast -m span --bp_bridge 50 --threads 4
```

Note that if you have fewer than 4 threads, you should change this accordingly. Increasing past 4 threads is not advantageous, as it primarily affects the BLAST step.

Three new directories have been written to the `/bin/05-Reference-blast/Output/` folder, which contains all the outputs.

The main results are the filtered fasta files, which are contained in the `/bin/05-Reference-blast/Output/Filtered-Fasta-Files/` directory. 

These filtered fasta files are ready for the sequence selection step.

---------------

### **4. Sequence Selection**

Here we select a single sequence per species per locus (Arthroleptids), and select a single sequence per sample per locus (Hyperoliids). We perform the second option by using the "subspecies" option that allows for voucher codes to be used as the third part of a three-part taxon name.

#### **Select Representative Sequences**

For supermatrices, a single sequence is typically used to represent each species for each gene. If multiple sequences for a given gene exist for a given species (e.g., because multiple individuals were sampled), then an objective strategy must be used for sequence selection. We can use the `Filter_Seqs_and_Species.py` module to accomplish this. 

We will select sequences for all loci based only on length (longest). We will also set 200 bp as the minimum length required to keep a sequence.

To prepare for this step, you will need to move or copy ALL of the similarity filtered fasta files from the previous step. There should be 6 filtered fasta files in total. This includes all the files in the `/bin/05-Reference-blast/Output/Filtered-Fasta-Files/` directory. Go ahead and move or copy these files to the `/bin/06-Filter-seqs/Input/` folder.

To run `Filter_Seqs_and_Species.py`, we need the filtered fasta files and a taxon names list file (`Taxon-List.txt`).

We will run `Filter_Seqs_and_Species.py` using the following command:

```
python Filter_Seqs_and_Species.py -i /bin/06-Filter-seqs/Input/ -o /bin/06-Filter-seqs/Output/ -s oneseq -f length -m 200 -t /bin/03-Starting-materials/Taxon-List.txt
```

The outputs are written to three folders contained in `/bin/06-Filter-seqs/Output/Results/`. 

The main outputs are the fasta files with one sequence per taxon, which are located in `/bin/06-Filter-seqs/Output/Results/Filtered-Fasta-Files/`. We will use these for the next step.

---------------

### **5. Sequence Alignment**

In preparation for alignment, we need to complete one pre-alignment step. We will adjust the direction of all sequences within each fasta file. This will ensure they are all oriented in the same direction, and prevent catastrophes in sequence alignment.

#### **Pre-Alignment: Adjust sequence directions**

We will ensure the sequences in each fasta file are correctly oriented by using the `Adjust_Direction.py` module. All of the target fasta files should be in the same output directory after the sequence selection step, so we will run this step using that location. 

Here is the command we will use:

```
python Adjust_Direction.py -i /bin/06-Filter-seqs/Output/Results/Filtered-Fasta-Files/ -o /bin/07-Adjust/ --threads 8
```

You should adjust the number of threads accordingly.

There should now be two output directories and a summary file in the `/bin/07-Adjust/` folder. The adjusted fasta files are the main output we are interested in, and they are located in the `/bin/07-Adjust/Adjusted-Fasta-Files/` folder. We will use these for alignment.

#### **Multiple Sequence Alignment**

We will use MAFFT (--auto) to align all 6 loci. This is by far the fastest alignment option, and should take ~1 minute total with 8 threads.

We will perform sequence alignments using the `Align.py` module using the following command:

```
python Align.py -i /bin/07-Adjust/Adjusted-Fasta-Files/ -o /bin/08-Align/ -a mafft --threads 8
```

You should adjust the number of threads accordingly.

The final alignments will be written to `/bin/08-Align/Alignments-MAFFT/`.

---------------

### **6. Post-Alignment**

We will now perform two post-alignment tasks:

- We will relabel the sequence records with relevant taxon names, making them compatible with downstream programs.
- We will concatenate the alignments to make a final supermatrix.

#### **Relabel sequences in alignments**

Currently, the alignments contain sequences that are labeled using a modified version of the GenBank record line. This cannot be used for most downstream programs. We will relabel these sequences to include the taxon label. In this case, we are treating the voucher codes for the hyperoliid sequences as a component of a three-part taxon name, so these will be written too if we use the subspecies option.

All of the target fasta files are in the `/bin/08-Align/Alignments-MAFFT/` folder.

We will relabel sequences with the `Fasta_Relabel_Seqs.py` module using the following command:

```
python Fasta_Relabel_Seqs.py -i /bin/08-Align/Alignments-MAFFT/ -o /bin/09-Relabel/ -r species -s /bin/03-Starting-materials/Taxon-List.txt
```

Here, the `-s` flag indicates we want to use the subspecies names present in the `Taxon-List.txt` file. These include the three-part names that incorporate a voucher code, so we can write this style of name for the hyperoliid sequences. The arthroleptid names consist only of two parts (genus species), so they will be written in the species style.

As a result, the sequence labels vary in this dataset because of this strategy of dealing with voucher codes (treating them as part of the taxon name). They include a mix of species + voucher and species labels, such as `Afrixalus_wittei_UTEP21795` and `Arthroleptis_adelphus`. If you open one of the relabeled alignments, this will become immediately apparent.

The fasta alignments with relabeled sequences will be written to `/bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/`.


#### **Concatenate alignments to create the final supermatrix**

The moment is upon us! It is time to create the final supermatrix. We will use the `Concatenation.py` module to perform this task. We will use the relabeled fasta files for this step, which are in the `/bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/` directory. 

We will specify the input format is fasta, the output format should be phylip (compatible with RAxML), and that missing characters should be coded with a `-` symbol.

We will perform this task using the following command:

```
python Concatenation.py -i /bin/09-Relabel/Relabeled-by-species/Relabeled-fasta-files/ -o /bin/10-Concatenate/ --informat fasta --outformat phylip -s dash
```

And here is the output on screen when we run it:

```
Found 6 fasta files to concatenate.

Found 365 unique taxa across alignment files.

Gathering sequences for all taxa (this could take some time)...
	Done.

Writing partitions file.

Writing concatenated alignment.
	Total alignment length = 4,550 bp.
	Total number of sequences included = 1,734.

Finished. Total elapsed time: 0:00:00.045934 (H:M:S)
```

A concatenated alignment in phylip format (`Concatenated_Alignment.phy`), a data partitions file (`Data_Partitions.txt`), and final count of the number of taxa and their number of corresponding loci (`Taxa_Loci_Count.log`) will be written to `/bin/10-Concatenate/`. 

And that's it! We've just successfully incorporated a large number of GenBank sequences into our local supermatrix project! I encourage you to look at the `Taxa_Loci_Count.log` to get a sense for how this worked. We added GenBank species sampling to our local vouchered population-level sampling to create a custom dataset suited for our project.

You can use the concatenated alignment to create a phylogeny using RAxML, like we did in our study. Our bootstrapped tree is available in the `OSF/11-RAxML/` folder, in case you'd like to download it. Again, you can see how hyperoliids are represented by voucher-level sampling, and the outgroup is represented by species-level sampling. That's pretty useful!


