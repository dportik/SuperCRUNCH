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

#### Basic Usage:

```
python Remove_Duplicate_Accessions.py -i <fasta file> -o <output directory>
```

#### Argument Explanations:

##### `-i <full-path-to-file>` 
 
> **Required**:  The full path to a fasta file with GenBank sequence data to filter.

##### `-o <path-to-directory>` 

> **Required**: The full path to an existing directory to write  output fasta file.





## **Taxon Filtering and Locus Extraction** <a name="TFLE"></a>

![F1](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig1.jpg)

Something




### Taxa_Assessment.py <a name="TA"></a>

Something

#### Basic Usage:

```
python Taxa_Assessment.py -i <fasta file> -t <taxon file> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> **Required**: The full path to a fasta file of GenBank sequence data.

##### `-t <path-to-file>`

> **Required**: The full path to a text file containing all taxon names  to cross-reference in the fasta file.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write  output files.

##### `--no_subspecies`

> **Optional**: Ignore subspecies labels in both the taxon names file and the fasta file.




### Rename_Merge.py <a name="RM"></a>

Something

#### Basic Usage:

```
python Rename_Merge.py -i <fasta file> -r <taxon renaming file> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> **Required**: The full path to a fasta file with taxon names to replace (*Unmatched_Taxa.fasta*).

##### `-r <path-to-file>`

> **Required**: The full path to a two-column text file containing all taxon names to be replaced, and the replacement names.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write output files.

##### `-m <full-path-to-file>`

> **Optional**: The full path to a fasta file containing valid taxon names (*Matched_Taxa.fasta*). 




### Parse_Loci.py <a name="PL"></a>

Something

#### Basic Usage:

```
python Parse_Loci.py -i <fasta file> -l <locus term file> -t <taxon file> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> **Required**: The full path to a fasta file of GenBank sequence data.

##### `-l <path-to-file>`

> **Required**: The full path to a three-column text file containing loci information to search for within the fasta file.

##### `-t <path-to-file>`

> **Required**: The full path to a text file containing all taxon names to cross-reference in the fasta file.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write output files.

##### `--no_subspecies`

> **Optional**: Ignore subspecies labels in both the taxon names file and the fasta file.






## **Orthology Filtering** <a name="OF"></a>

![F2](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig2.jpg)

Something




### Cluster_Blast_Extract.py <a name="CBE"></a>

Something

#### Basic Usage:

```
python Cluster_Blast_Extract.py -i <fasta file directory> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the parsed locus-specific fasta files.

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. Recommended: *dc-megablast*.

##### `-m <choice>`

> **Optional**: The strategy for dealing with multiple non-overlapping blast coordinates. Choices = *span, nospan, all*. Default = *span*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets.




### Reference_Blast_Extract.py <a name="RBE"></a>

Something

#### Basic Usage:

```
python Reference_Blast_Extract.py -i <input directory> -d <reference fasta name> -e <empirical fasta name> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the reference fasta file and the empirical fasta file.

##### `-d <filename>`

> **Required**: The name of the reference fasta file that will be used to create the blast database. Requires file name only, NOT full path, as it should be located in the input directory (-i).

##### `-e <filename>`

> **Required**: The name of the empirical fasta file to blast to the database to prune sequences. Requires file name only, NOT full path, as it should be located in the input directory (-i).

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. Recommended: *dc-megablast*.

##### `-m <choice>`

> **Optional**: The strategy for dealing with multiple non-overlapping blast coordinates. Choices = *span, nospan, all*. Default = *span*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets.




### Contamination_Filter.py <a name="CF"></a>

Something

#### Basic Usage:

```
python Contamination_Filter.py -i <input directory> -d <contamination fasta name> -e <empirical fasta name> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the reference fasta file and the empirical fasta file.

##### `-d <filename>`

> **Required**: The name of the contamination reference fasta file that will be used to create the blast database. Requires file name only, NOT full path, as it should be located in the input directory (-i).

##### `-e <filename>`

> **Required**: The name of the empirical fasta file to blast to the database to find bad sequences. Requires file name only, NOT full path, as it should be located in the input directory (-i).

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. Recommended: *megablast*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence.






## **Sequence Quality Filtering and Selection** <a name="SQFS"></a>

![F3](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig3.jpg)

Something




### Filter_Seqs_and_Species.py <a name="FSS"></a>

Something

#### Basic Usage:

```
python Filter_Seqs_and_Species.py -i <input directory> -f <filter strategy> -l <minimum base pairs> -t <taxon file>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter.

##### `-f <choice>`

> **Required**: Strategy for filtering sequence data. Choices = *translate, length*.

##### `-l <integer>`

> **Required**: An integer for the minimum number of base pairs required to keep a sequence (ex. 150).

##### `-t <path-to-file>`

> **Required**: The full path to a text file containing all taxon names to cross-reference in the fasta file.

##### `--table <choice>`

> **Required for** `-f translate`: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*.

##### `--randomize`

> **Optional**: For taxa with multiple sequences, shuffle order randomly. Overrides sorting by length for all methods (-f).

##### `--allseqs`

> **Optional**: For taxa with multiple sequences, select all sequences passing the filters instead of a single representative sequence.

##### `--no_subspecies`

> **Optional**: Ignore subspecies labels in both the taxon names file and the fasta file.




### Make_Acc_Table.py <a name="MAT"></a>

Something

#### Basic Usage:

```
python Make_Acc_Table.py -i <input directory>
```

#### Argument Explanations:

##### -i <path-to-directory>

> **Required**: The full path to a directory containing the fasta files (single sequence per taxon).

##### -s <path-to-file>

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file.







## **Sequence Alignment** <a name="SA"></a>

![F4](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig4.jpg)

Something




### Adjust_Direction.py <a name="AD"></a>

Something

#### Basic Usage:

```
python Adjust_Direction.py -i <input directory>
```

#### Argument Explanations:

##### -i <path-to-directory>

> **Required**: The full path to a directory which contains the unaligned fasta files.




### Coding_Translation_Tests.py <a name="CTT"></a>

Something

#### Basic Usage:

```
python Coding_Translation_Tests.py -i <input directory> --table <translation table>
```

##### -i <path-to-directory> 

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter.

##### --table <choice>

> **Required**: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*.

##### --rc

> **Optional**: In addition to forward frames, use reverse complement for translation tests. Not recommended if direction of sequences has already been adjusted.




### Align.py <a name="A"></a>

Something

#### Basic Usage:

```
python Align.py -i <input directory> -a <aligner> 
```

##### -i <path-to-directory>

> **Required**: The full path to a directory which contains the unaligned fasta files.

##### -a <choice>

> **Required**: Specify whether alignment is by mafft, macse, muscle, or clustalo. If macse, MUST provide flags --mpath and --table. Choices = *mafft, macse, muscle, clustalo, all*.

##### --mpath <path-to-executable>

> **Required** for `-a macse`: Full path to a macse jar file.

##### --table <choice>

> **Required** for `-a macse`: Specifies translation table. Choices = *standard, vmtdna*.

##### --mem <integer>

> **Optional** for `-a macse`: An integer for how much memory to assign to macse (in GB). Default = 1.

##### --pass_fail

> **Optional** for `-a macse`: Specifies macse to perform dual alignment. Files in -i directory must follow labeling format: NAME_Passed.fasta, NAME_Failed.fasta.

##### --accurate

> **Optional**: Specifies more thorough search settings in mafft, clustalo, or macse (see below).




### Trim_Alignments.py <a name="TAS"></a>

Something

#### Basic Usage:

```
python Trim_Alignments.py -i <input directory> -f <output format> -a <trimal method>
```

##### -i <path-to-directory>

> **Required**: The full path to a directory which contains the input alignment files. File formats can include fasta, nexus, or phylip (see below).

##### -f <choice>

> **Required**: Specifies the output file format for trimmed alignments. Choices = *fasta, nexus, phylip*.

##### -a <choice>

> **Required**: Specifies the trimal method for trimming alignments. Choices = *gt, noallgaps, both*.

##### --gt <value>

> **Optional**: Specifies the gap threshold (gt) value for trimal, the minimum fraction of sequences without a gap. Must be between 0 and 1. Default = 0.05.







## **File Formatting Tasks** <a name="FFT"></a>

![F5](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig5.jpg)

Something




### Relabel_Fasta.py <a name="RF"></a>

Something

#### Basic Usage:

```
python Relabel_Fasta.py -i <input directory> -r <relabel option>
```

##### -i <path-to-directory>

> **Required**: The full path to a directory containing the unaligned or aligned fasta files.

##### -r <choice>

> **Required**: The strategy for relabeling sequence records. Choices = *species, accession, species_acc*.

##### -s <path-to-file>

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file.




### Fasta_Convert.py <a name="FC"></a>

Something

#### Basic Usage:

```
python Fasta_Convert.py -i <input directory>
```

##### -i <path-to-directory>

> **Required**: The full path to a directory which contains the aligned fasta files.




### Concatenation.py <a name="C"></a>

Something

#### Basic Usage:

```
python Concatenation.py -i <input directory> -r <input format> -s <missing data symbol> -o <output format>
```

##### -i <path-to-directory>

> **Required**: The full path to a directory containing the aligned files.

##### -f <choice>

> **Required**: The input file format of the alignments. Choices = *fasta, phylip*.

##### -s <choice>

> **Required**: A base pair symbol used to represent missing data sequences. Choices = *dash, N, ?*.

##### -o <choice>

> **Required**: The output file format for the final concatenated alignment. Choices = *fasta, phylip*.



The end.

------
*Last updated: January 2019*
