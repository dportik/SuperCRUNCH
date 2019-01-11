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

#### Arguments:
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

### Rename_Merge.py <a name="RM"></a>

Something

### Parse_Loci.py <a name="PL"></a>

Something

## **Orthology Filtering** <a name="OF"></a>

![F2](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig2.jpg)

Something

### Cluster_Blast_Extract.py <a name="CBE"></a>

Something

### Reference_Blast_Extract.py <a name="RBE"></a>

Something

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
