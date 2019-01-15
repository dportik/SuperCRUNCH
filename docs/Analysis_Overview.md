![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Analysis Overview

Complete instructions for performing each step of **SuperCRUNCH** are provided here. The scripts are listed in the approximate order in which they should be used, and split among larger topics of interest. Helpful information can be displayed for all scripts on the command line by running them using only the -h flag. 

### [Starting Materials](#GSM):
+ [Obtaining Sequence Data](#OSD)
+ Optional tool: [Remove_Duplicate_Accessions](#RDA)
+ [Obtaining Taxon Names Lists](#OTNL)
+ [Obtaining Loci and Search Terms](#OLST)

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




---------------

## **Starting Material** <a name="GSM"></a>

To run **SuperCRUNCH**, you will need a fasta file of downloaded nucleotide sequence records, a list of taxa, and a list of loci with associated search terms.

### Obtaining Sequence Data <a name="OSD"></a>

The starting molecular data for **SuperCRUNCH** consists of a fasta file of downloaded nucleotide sequence records from the [NCBI nucleotide database (GenBank)](https://www.ncbi.nlm.nih.gov/nucleotide/). The simplest way to obtain these data is to search for a taxonomic term on the nucleotide database and download all the available records in fasta format. However, for larger taxonomic groups (ex. all birds or frogs) this may not be possible as there will be too much data to download directly. In this case, it is better to split searches using the taxonomy of the group (order, clade, family, etc), download each record set in fasta format, and then combine the resulting fasta files. The [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) is particularly useful for identifying the relevant taxonomic search terms. Advanced users may choose to use additional search terms to filter the results, which could help reduce the resulting file size, but this is not necessary. The downloaded fasta files may be quite large in size (several GB) and it is unlikely you can use a text editor to open and edit them, so combining files requires a different approach. The fasta files will all have the extension '.fasta', and they can be quickly combined using Unix. If the fasta files are all in the same working directory, you can use the following command to accomplish this:

`cat *.fasta > Combined.fasta`

The starting fasta file can be in the interleaved format or sequential format, as demonstrated below: 

**Interleaved format:**
```
>KP820543.1 Callisaurus draconoides voucher MVZ265543 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATTTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTGGGTAATGTAATGGCGAGCGAT
TGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATTCTCAAACAAGAGCAGACAGAGTGCCCCCAGG

>KP820544.1 Urosaurus ornatus voucher UWBM7587 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATCTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTAGGTAATGTAATGGCCAGCGAT
TGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATCCTCAAACAAGAACAGACTGAGTGCCCCCAGG
```

**Sequential format:**
```
>KP820543.1 Callisaurus draconoides voucher MVZ265543 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATTTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTGGGTAATGTAATGGCGAGCGATTGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATTCTCAAACAAGAGCAGACAGAGTGCCCCCAGG

>KP820544.1 Urosaurus ornatus voucher UWBM7587 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATCTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTAGGTAATGTAATGGCCAGCGATTGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATCCTCAAACAAGAACAGACTGAGTGCCCCCAGG
```

The record labels must have a > with a unique accession number, followed by a sequence description. If curated properly, the description line will contain a taxon name, locus information, and possibly several other identifiers, all separated by whitespace. There should not be any vertical bars (|) separating information contained in the sequence label. **SuperCRUNCH** cannot properly parse this format.

The following vertical bar fasta format is ***NOT*** supported by **SuperCRUNCH**:

```
>gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA
CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAGAATATATGATCGAGTG
AATCTGGAGGACCTGTGGTAACTCAGCTCGTCGTGGCACTGCTTTTGTCGTGACCCTGCTTTGTTGTTGG

>gi|2765656|emb|Z78531.1|CFZ78531 C.fasciculatum 5.8S rRNA gene and ITS1 and ITS2 DNA
CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAGAACATACGATCGAGTG
AATCCGGAGGACCCGTGGTTACACGGCTCACCGTGGCTTTGCTCTCGTGGTGAACCCGGTTTGCGACCGG
```

If your starting sequence set is in this format, you will need to convert it to the simpler whitespace fasta format before attempting to run **SuperCRUNCH**.

### Remove_Duplicate_Accessions.py <a name="RDA"></a>

Sometimes combining several fasta files can produce duplicate entries, in other words records with identical accession numbers. This will cause problems with the way **SuperCRUNCH** loads fasta files, as there can be no duplicate accession numbers. You will know this is the case if you try to run a script and it produces the following error:

`ValueError: Duplicate key 'AF443291.2'`

In this instance, there are two records each containing the accession number AF443291.2. In order to use **SuperCRUNCH**, you'll need to clean the fasta file to remove any duplicate accession numbers. This script can be used for this purpose.


#### Basic Usage:

```
python Remove_Duplicate_Accessions.py -i <fasta file> -o <output directory>
```

#### Argument Explanations:

##### `-i <full-path-to-file>` 
 
> **Required**:  The full path to a fasta file with GenBank sequence data to filter.

##### `-o <path-to-directory>` 

> **Required**: The full path to an existing directory to write  output fasta file.

By running this script, all duplicate entries will be removed and a new fasta file will be written to the output directory called *Cleaned.fasta*. This fasta file should be used to run **SuperCRUNCH**.

### Obtaining Taxon Names Lists <a name="OTNL"></a>

**SuperCRUNCH** requires a list of taxon names that is used to filter sequences. Lists of taxon names can be obtained through general databases, such as the NCBI Taxonomy Browser. In many cases there are specific databases dedicated to major groups, for example the [Reptile Database](http://www.reptile-database.org/), [AmphibiaWeb](https://amphibiaweb.org/), and [Amphibian Species of the World](http://research.amnh.org/vz/herpetology/amphibia/), which usually contain up-to-date taxonomies in a downloadable format. 

The taxon names list required is a simple text file which contains one taxon name per line. The file can contain a mix of species (binomial) and subspecies (trinomial) names, and components of each taxon name should be separated by a space (rather than undescore). Sometimes using Excel or other applications to generate the taxon names text file will include hidden characters not compatible with **SuperCRUNCH**, so make sure to open the file in a text editor and ensure the format includes Unix line breaks and UTF-8 encoding.
Below are some example contents from suitable taxon name lists:

A partial list containing species (binomial) names:

```
Lygodactylus regulus
Lygodactylus rex
Lygodactylus roavolana
Lygodactylus scheffleri
Lygodactylus scorteccii
Lygodactylus somalicus
```

A partial list containing species (binomial) and subspecies (trinomial) names:

```
Leycesteria crocothyrsos
Leycesteria formosa
Linnaea borealis americana
Linnaea borealis borealis
Linnaea borealis longiflora
```

**SuperCRUNCH** offers the option to exclude subspecies from searches, so a taxon list containing a mix of species and subspecies does not need to be pruned by hand if subspecies are not desired in the analysis. The effect of this option with different types of taxon lists is explained in usage of the *Taxa_Assessment.py* script.

### Obtaining Loci and Search Terms <a name="OLST"></a>

**SuperCRUNCH** requires a list of loci and associated search terms to initially identify the content of sequence records. For each locus included in the list, **SuperCRUNCH** will search for an associated abbreviated name and full description in the record label. The choice of loci to include will inherently be group-specific, and surveys of phylogenetic and phylogeographic papers may help identify an appropriate marker set. There is no limit to the number of loci, and **SuperCRUNCH** can also be used to search for large genomic data sets, such as those obtained through sequence capture experiments (UCEs, anchored enrichment, etc.). 

The format of the locus text file involves three tab-delimited columns. The first column should contain the locus name that will be used to label output files. It should not contain any spaces or special characters. The second column should contain the known abbreviation for the gene or marker. This second column can contain multiple abbreviations, which should be separated by a semi-colon (with no extra spaces between names). The third column should contain the full label of the gene or marker. This third column can also contain multiple search entries, which should be separated by a semi-colon (with no extra spaces between entries). The abbreviations and labels are not case-specific, and they are all converted to uppercase during actual searches (along with the sequence record labels).

Here is an example of the formatting for a locus file containing three genes to search for:

```
CMOS	CMOS;C-MOS	oocyte maturation factor
EXPH5	EXPH5	exophilin;exophilin 5;exophilin-5;exophilin protein 5
PTPN	PTPN;PTPN12	protein tyrosine phosphatase;tyrosine phosphatase non-receptor type 12
```

In this file:
+ CMOS contains two abbreviations and one label search term. 
+ EXPH5 contains one abbreviation and four label search terms.
+ PTPN contains two abbreviations and two label search terms. 

How does the actual locus searching work? 

For locus abbreviations, the sequence record label is split by spaces, stripped of punctuation, and converted to uppercase. Each resulting component is checked to see if it is identical to an included locus abbreviation. If so, a match is found. 

> If the locus file contains:

`CMOS	cmos;c-mos	oocyte maturation factor`

> The abbreviations will include:

```
CMOS
C-MOS
```

> Notice the search terms are converted to uppercase

> If the sequence record contains the following label:

`>JX838886.1 Acanthocercus annectens voucher CAS 227508; oocyte maturation factor (CMOS) gene, partial cds`

> Then it will result in the following components:

```
>JX838886.1
ACANTHOCERCUS
ANNECTENS
VOUCHER
CAS
227508
OOCYTE
MATURATION
FACTOR
CMOS
GENE
PARTIAL
CDS
```

> Notice the line is stripped of all punctuation (including parentheses) and converted to uppercase. In this example, a match will be found using the CMOS search term, but not the C-MOS term. 

For locus labels, the sequence record label is converted to uppercase (punctuation is left intact). The line is then checked to see if contains an included locus label. If so, a match is found. 

> If the locus file contains:

`EXPH5	EXPH5	exophilin;exophilin 5;exophilin-5;exophilin protein 5`

> The labels will include:

```
EXOPHILIN
EXOPHILIN 5
EXOPHILIN-5
EXOPHILIN PROTEIN 5
```

> If the sequence record contains the following label:

`>JX999516.1 Liolaemus pictus voucher LP111; exophilin 5 (EXPH5) gene, partial cds`

> It will be converted to the following search line:

`>JX999516.1 LIOLAEMUS PICTUS VOUCHER LP111; EXOPHILIN 5 (EXPH5) GENE, PARTIAL CDS`

> Notice punctuation and parentheses are left intact. In this example, the EXOPHILIN 5 and EXOPHILIN labels will both find a match, but not EXOPHILIN-5 or EXOPHILIN PROTEIN 5.

The success of finding loci depends on defining appropriate locus abbreviations and labels. 

---------------

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





---------------

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




---------------

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





---------------

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





---------------

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
