![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Analysis Overview

Complete instructions for performing each step of **SuperCRUNCH** are provided here. The scripts are listed in the approximate order in which they should be used, and split among larger topics of interest.

### [Getting the Starting Material](#GSM):
+ Optional tool: [Remove_Duplicate_Accessions](#RDA)

### [Taxon Filtering and Locus Extraction](#TFLE):

+ [Taxa_Assessment](#TA)
+ [Rename_Merge](#RM)
+ [Parse_Loci](#PL)

### [Orthology Filtering](#OF)

+ [Cluster_Blast_Extract](#CBE)
+ [Reference_Blast_Extract](#RBE)
+ [Contamination_Filter](#CF)

### [Sequence Quality Filtering and Selection](#SQFS)

+ [Filter_Seqs_and_Species](#FSS)
+ [Make_Acc_Table](#MAT)

### [Sequence Alignment](#SA)

+ [Adjust_Direction](#AD)
+ [Coding_Translation_Tests](#CTT)
+ [Align](#A)
+ [Trim_Alignments](#TAS)

### [File Formatting Tasks](#FFT)

+ [Relabel_Fasta](#RF)
+ [Fasta_Convert](#FC)
+ [Concatenation](#C)

The steps involved in a typical run include executing a majority of these steps. However, there is a lot of flexibility and  epending on the analysis goals.   

## **Getting the Starting Material** <a name="GSM"></a>

Something

### Remove_Duplicate_Accessions.py <a name="RDA"></a>

Something.

#### Usage:
| Argument | Type | Description |
| ----- | ----- | ----- | 
| -i | Required | The full path to a fasta file with GenBank sequence data to filter. |
| -o | Required | The full path to an existing directory to write output fasta file. | 
> ``
> ``


## **Taxon Filtering and Locus Extraction** <a name="TFLE"></a>

![F1](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig1.jpg)

Something


### Taxa_Assessment.py <a name="TA"></a>

Something

#### Usage:
| Argument | Type | Description |
| ---------------- | ----- | ----- | 
| -i | Required | The full path to a fasta file of GenBank sequence data. |
| -t | Required | The full path to a text file containing all taxon names to cross-reference in the fasta file. | 
|  -o  | Required | The full path to an existing directory to write output files. | 
|  -- no_subspecies  | *Optional* | Ignore subspecies labels in both the taxon names file and the fasta file. | 

### Rename_Merge.py <a name="RM"></a>

Something

#### Usage:
| Argument | Type | Description |
| ----- | ----- | ----- | 
| -i | Required | The full path to a fasta file with taxon names to replace (*Unmatched_Taxa.fasta*). |
| -r | Required | The full path to a two-column text file containing all taxon names to be replaced, and the replacement names. | 
| -o | Required | The full path to an existing directory to write output files. | 
| -m | *Optional* | The full path to a fasta file containing valid taxon names (*Matched_Taxa.fasta*). | 


### Parse_Loci.py <a name="PL"></a>

Something

#### Usage:
| Argument | Type | Description |
| ---------------- | ----- | ----- | 
| -i | Required | The full path to a fasta file of GenBank sequence data. |
| -l | Required | The full path to a three-column text file containing loci information to search for within the fasta file. | 
| -t | Required | The full path to a text file containing all taxon names to cross-reference in the fasta file. | 
| -o | Required | The full path to an existing directory to write output files. | 
| --no_subspecies | *Optional* | Ignore subspecies labels in both the taxon names file and the fasta file. | 

## **Orthology Filtering** <a name="OF"></a>

![F2](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig2.jpg)

Something

### Cluster_Blast_Extract.py <a name="CBE"></a>

Something

#### Usage:
| Argument | Type | Description |
| ---------- | ----- | ----- | 
| -i | Required | The full path to a directory containing the parsed locus-specific fasta files. |
| -b | Required | The blast algorithm to use. Options = blastn, blastn-short, dc-megablast, megablast. | 
| -m | *Optional* | The strategy for dealing with multiple non-overlapping blast coordinates. Options = span, nospan, all. Default = span. | 
| --max_hits | *Optional* | The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets. | 

### Reference_Blast_Extract.py <a name="RBE"></a>

Something

#### Usage:
| Argument | Type | Description |
| ---------- | ----- | ----- | 
| -i | Required | The full path to a directory containing the reference fasta file and the empirical fasta file. |
| -d | Required | The name of the reference fasta file that will be used to create the blast database. Requires file name only, NOT full path, as it should be located in the input directory (-i). |
| -e | Required | The name of the empirical fasta file to blast to the database to prune sequences. Requires file name only, NOT full path, as it should be located in the input directory (-i). | 
| -b | Required | The blast algorithm to use. Options = blastn, blastn-short, dc-megablast, megablast. | 
| -m | *Optional* | The strategy for dealing with multiple non-overlapping blast coordinates. Options = span, nospan, all. Default = span. | 
| --max_hits | *Optional* | The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets. | 

### Contamination_Filter.py <a name="CF"></a>

Something

## **Sequence Quality Filtering and Selection** <a name="SQFS"></a>

![F3](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig3.jpg)

Something

### Filter_Seqs_and_Species.py <a name="FSS"></a>

Something

### Make_Acc_Table.py <a name="MAT"></a>

Something

## **Sequence Alignment** <a name="SA"></a>

![F4](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig4.jpg)

Something

### Adjust_Direction.py <a name="AD"></a>

Something

### Coding_Translation_Tests.py <a name="CTT"></a>

Something

### Align.py <a name="A"></a>

Something

### Trim_Alignments.py <a name="TAS"></a>

Something

## **File Formatting Tasks** <a name="FFT"></a>

![F5](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig5.jpg)

Something

### Relabel_Fasta.py <a name="RF"></a>

Something

### Fasta_Convert.py <a name="FC"></a>

Something

### Concatenation.py <a name="C"></a>

Something
