![SuperCrunch Logo](https://github.com/dportik/SuperCRUNCH/blob/master/docs/SuperCRUNCH_Logo.png)

---------------

## Analysis Overview and Detailed Instructions

Complete instructions for performing each step of **SuperCRUNCH** are provided here. The scripts are listed in the approximate order in which they should be used, and split among larger topics of interest. Helpful information can also be displayed for all scripts on the command line by running them using the -h flag. 

### [Starting Materials](#GSM):
+ [Obtaining Sequence Data](#OSD)
+ Optional tool: [Remove_Duplicate_Accessions](#RDA)
+ [Obtaining Taxon Names Lists](#OTNL)
+ [Obtaining Loci and Search Terms](#OLST)
    + [Searching for UCE loci](#SFUL)

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

The steps involved in a typical run include executing a majority of these steps. However, there is a lot of flexibility and the workflow can be tailored to achieve various analysis goals.   




---------------

## **Starting Material** <a name="GSM"></a>

To run **SuperCRUNCH**, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms. Information on how to create these files is below.

---------------

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

---------------

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

---------------

### Obtaining Taxon Names Lists <a name="OTNL"></a>

**SuperCRUNCH** requires a list of taxon names that is used to filter sequences. Lists of taxon names can be obtained through general databases, such as the NCBI Taxonomy Browser. In many cases there are specific databases dedicated to major groups, for example the [Reptile Database](http://www.reptile-database.org/), [AmphibiaWeb](https://amphibiaweb.org/), and [Amphibian Species of the World](http://research.amnh.org/vz/herpetology/amphibia/), which usually contain up-to-date taxonomies in a downloadable format. 

The taxon names list required is a simple text file which contains one taxon name per line. The file can contain a mix of species (binomial) and subspecies (trinomial) names, and components of each taxon name should be separated by a space (rather than undescore). Sometimes using Excel or other applications to generate the taxon names text file will include hidden characters not compatible with **SuperCRUNCH**, so make sure to open the file in a text editor and ensure the format includes Unix line breaks (line breaks marked by \n, rather than \r\n) and UTF-8 encoding. Below are some example contents from suitable taxon name lists.

Partial contents of a list containing species (binomial) names:

```
Lygodactylus regulus
Lygodactylus rex
Lygodactylus roavolana
Lygodactylus scheffleri
Lygodactylus scorteccii
Lygodactylus somalicus
```

Partial contents of a list containing species (binomial) and subspecies (trinomial) names:

```
Leycesteria crocothyrsos
Leycesteria formosa
Linnaea borealis americana
Linnaea borealis borealis
Linnaea borealis longiflora
```

**SuperCRUNCH** offers the option to exclude subspecies from searches, so a taxon list containing a mix of species and subspecies does not need to be pruned by hand if subspecies are not desired in the analysis. The effects of this option with different types of taxon lists is explained in usage of the `Taxa_Assessment.py` script.

---------------

### Obtaining Loci and Search Terms <a name="OLST"></a>

**SuperCRUNCH** requires a list of loci and associated search terms to initially identify the content of sequence records. For each locus included in the list, **SuperCRUNCH** will search for an associated abbreviated name and longer label in the sequence record. The choice of loci to include will inherently be group-specific, and surveys of phylogenetic and phylogeographic papers may help to identify an appropriate marker set. There is no limit to the number of loci, and **SuperCRUNCH** can also be used to search for large genomic data sets, such as those obtained through sequence capture experiments (UCEs, anchored enrichment, etc.). For detailed instructions on searching for UCEs, see the next section.

The format of the locus text file involves three tab-delimited columns. The first column contains the locus name that will be used to label output files. It must not contain any spaces or special characters (including underscores and hyphens), and should be kept short and simple. The second column contains the known abbreviation(s) for the gene or marker. Abbreviations should not include any spaces, but can contain special characters. The second column can contain multiple abbreviations, which should be separated by a semi-colon (with no extra spaces between names). The third column contains a longer label of the gene or marker, such as its full name or description. This third column can also contain multiple search entries, which should be separated by a semi-colon (with no extra spaces between entries). The abbreviations and labels are not case-specific, and they are all converted to uppercase during actual searches (along with the sequence record labels). Although the locus text file can be created using Excel or other applications, it must be a tab-delimited text file. Similar to the taxon names file, make sure to open the file in a text editor and ensure the format includes Unix line breaks (line breaks marked by \n, rather than \r\n) and UTF-8 encoding, otherwise extra characters may interfere with parsing the file correctly with **SuperCRUNCH**.

The success of finding loci depends on defining appropriate locus abbreviations and labels. Examples of how searches work can be found below, which should guide how to select good search terms. For any locus, it is a good idea to examine several records using GenBank to identify commonly used labels.

Here is an example of the formatting for a locus file containing three genes to search for:

```
CMOS	CMOS;C-MOS	oocyte maturation factor
EXPH5	EXPH5	exophilin;exophilin 5;exophilin-5;exophilin protein 5
PTPN	PTPN;PTPN12	protein tyrosine phosphatase;tyrosine phosphatase non-receptor type 12
```

These are the partial contents of a locus search terms file I used for squamates (reptiles). The full file (*Locus-Search-Terms_Squamate_Markers.txt*) is provided in the example data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data).

In this example:
+ CMOS contains two abbreviations and one label search term. 
+ EXPH5 contains one abbreviation and four label search terms.
+ PTPN contains two abbreviations and two label search terms. 

**How does the actual locus searching work?**

+ For locus abbreviations, the sequence record label is split by spaces, stripped of punctuation, and converted to uppercase. Each resulting component is checked to see if it is identical to an included locus abbreviation. If so, a match is found. 

+ For locus labels, the sequence record label is converted to uppercase (punctuation is left intact). The line is then checked to see if contains an included locus label. If so, a match is found. 

+ If a locus abbreviation *or* a locus label is matched to the contents of a sequence record, the record will pass the filtering step. 

**Example of locus abbreviation search:**

If the locus file contains:

`CMOS	cmos;c-mos	oocyte maturation factor`

The abbreviations will include:

```
CMOS
C-MOS
```

Notice the search terms are converted to uppercase

If the sequence record contains the following label:

`>JX838886.1 Acanthocercus annectens voucher CAS 227508; oocyte maturation factor (CMOS) gene, partial cds`

Then it will result in the following components:

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

Notice the line is stripped of all punctuation (including parentheses) and converted to uppercase. In this example, a match will be found using the CMOS search term, but not the C-MOS term.

**Example of locus label search:**

If the locus file contains:

`EXPH5	EXPH5	exophilin;exophilin 5;exophilin-5;exophilin protein 5`

The labels will include:

```
EXOPHILIN
EXOPHILIN 5
EXOPHILIN-5
EXOPHILIN PROTEIN 5
```

If the sequence record contains the following label:

`>JX999516.1 Liolaemus pictus voucher LP111; exophilin 5 (EXPH5) gene, partial cds`

It will be converted to the following search line:

`>JX999516.1 LIOLAEMUS PICTUS VOUCHER LP111; EXOPHILIN 5 (EXPH5) GENE, PARTIAL CDS`

Notice punctuation and parentheses are left intact and the line is simply converted to uppercase. The line is then checked to see if any supplied locus label is contained within. In this example, the `EXOPHILIN 5` and `EXOPHILIN` labels are both contained in the line and would produce a match, but `EXOPHILIN-5` and `EXOPHILIN PROTEIN 5` would not. The more specific or complex a label search term is, the less likely it is to produce an exact match. My recommendation is to find the simplest common denominator among records and include that label, along with more complex search labels.

---------------

### Searching for UCE loci <a name="SFUL"></a> 

The strategy for obtaining sets of UCE sequences is a little different from the smaller locus sets. First, you'll want to do a specific GenBank search for the taxonomic term of interest *and* `ultra conserved element` and/or `uce`. This will produce a much more manageable set of sequences to work with, as searching for several thousand loci in a large sequence set will inevitably take a very long time. 

To generate a locus search terms file, I retrieved the uce names from the ***uce-5k-probes.fasta*** file located [here](https://github.com/faircloth-lab/uce-probe-sets/tree/master/uce-5k-probe-set). Unfortunately, there does not appear to be a standard naming convention for the UCE loci on GenBank. If the sequences have been properly curated, then the description lines *should* contain the uce name somewhere (uce-10, uce-453, uce-5810, etc). If so, they will be compatible with the the 5k UCE locus search terms file (*Locus-Search-Terms_UCE_5k_set.txt*) I've made available in the data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data).

Here are partial contents from the UCE locus search term file:

```
uce-5805	uce-5805	xxxxxxxxxxxx
uce-5806	uce-5806	xxxxxxxxxxxx
uce-5808	uce-5808	xxxxxxxxxxxx
uce-5810	uce-5810	xxxxxxxxxxxx
```

Notice the third column is junk. Unfortunately, UCE loci have been numbered in a suboptimal way. For example, uce-1 instead of uce-0001. This causes problems when searching for locus labels, because the term `uce-1` is contained inside of `uce-10`, `uce-104`, `uce-1038`, etc. Because of this, the label search will not work properly, and so we have to rely exclusively on the abbreviation to find the correct sequences.

Here is an example of some frog UCE records I searched:

```
>KY160876.1 Kaloula kalingensis voucher RMB1887 ultra conserved element locus uce-5806 genomic sequence
ATATTTGTGTTTATTTTCTACTTGTATTAATTGACAACATTTGCCTGTTGGCTCAAGGGAATCAGTGTTC
CCATTTTATGCACTCTATTTTAAAATGCAGACAGTGGTAGAACAGATGTGTTTTTTTTAACCCCATA...

>KY160875.1 Kaloula pulchra voucher KU328278 ultra conserved element locus uce-5806 genomic sequence
ATATTTGTGTTTATTTTCTACTTGTATTAATTGACAACATTTGCCTGTTGGCTTAAGGGAATCATTGTTG
CCATTTTATGCACTCTATTTTAAAATGCATACAGTGGTAGAACAGATGTGTTTTTTTAACCCCATAG...

>KY160874.1 Kaloula picta voucher KU321376 ultra conserved element locus uce-5806 genomic sequence
ATATTTGTGTTTATTTTCTACTTGTATTAATTGACAACATTTGCCTGTTGGCTCAAGGGAATCAGTGTTG
CCATTTTATGCACTCTATTTTAAAATGCAGACAGTGGTAGAACAGATGTGTTTTTTTTAACCCCATA...

>KY160873.1 Kaloula conjuncta negrosensis voucher KU328639 ultra conserved element locus uce-5806 genomic sequence
ATATTTGTGTTTATTTTCTACTTGTATTAATTGACAACATTTGCCTGTTGGCTCAAGGGAATCAGTGTTG
CCATTTTATGCACTCTATTTTAAAATGCAGACAGTGGTAGAACAGATGTGTTTTTTTTAACCCCATA...
```

The locus abbreviation terms successfully retrieved the corresponding loci in the example above. 

The 5k UCE locus search terms file (*Locus-Search-Terms_UCE_5k_set.txt*) is freely available in the data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data), and it can be used to retrieve UCE data as long as records have the UCE locus name in the description lines.

---------------

## **Taxon Filtering and Locus Extraction** <a name="TFLE"></a>

![F1](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig1.jpg)

This section deals with the first filtering steps, which include screening sequence records for taxon names and gene/marker/locus identity. To run these steps in **SuperCRUNCH**, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms. Information on obtaining these files is provided in the section above. The `Taxa_Assessment.py` and `Rename_Merge.py` scripts are optional, but highly recommended. `Taxa_Assessment.py` identifies all valid and invalid taxon names contained within the starting fasta file. `Rename_Merge.py` is an optional data cleaning step that can be used to replace invalid taxon names with updated valid names for corresponding  sequence records. This relabeling step allows these records to pass the taxonomy filter in `Parse_Loci.py`, rather than be discarded. The original or updated fasta file is processed using `Parse_Loci.py`. For each locus included, a fasta file is produced. These fasta files contain records with valid taxon names that matched one or more of the search terms for the locus.

---------------

### Taxa_Assessment.py <a name="TA"></a>

The goal of this script is to search through records in a fasta file of NCBI nucleotide sequences to determine whether or not they contain a taxon name present in the user-supplied taxon list. The taxon names list can contain a mix of species (binomial name) and subspecies (trinomial name) labels, and searches are not case-sensitive. 

Two output fasta files are written to the specified output directory: one containing only records with valid taxon names (*Matched_Taxa.fasta*), and one containing records with invalid taxon names (*Unmatched_Taxa.fasta*). The accession numbers for each of these fasta files are also written to separate files (*Matched_Records_Accession_Numbers.log*, *Unmatched_Records_Accession_Numbers.log*). Two log files are written which contain lists of the valid (*Matched_Taxon_Names.log*) and invalid taxon names (*Unmatched_Taxon_Names.log*) found across all records. The *Unmatched_Taxon_Names.log* file can be used to create the file needed to relabel taxa in the `Rename_Merge.py` script. 

The decision to include or exclude subspecies labels is up to the user, and can be specified using the `--no_subspecies` flag. For a thorough explanation of how taxonomy searches are conducted and how this flag affects this step (and others), please see below.

#### Basic Usage:

```
python Taxa_Assessment.py -i <fasta file> -t <taxon file> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> **Required**: The full path to a fasta file of GenBank sequence data.

##### `-t <path-to-file>`

> **Required**: The full path to a text file containing all taxon names to cross-reference in the fasta file.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write  output files.

##### `--no_subspecies`

> **Optional**: Ignore the subspecies component of taxon labels in the taxon names file and in the fasta file to search.

#### Example Use:

```
python Taxa_Assessment.py -i bin/Analysis/Start_Seqs.fasta -t bin/Analysis/Taxa_List.txt -o bin/Analysis/Output/
```
```
python Taxa_Assessment.py -i bin/Analysis/Start_Seqs.fasta -t bin/Analysis/Taxa_List.txt -o bin/Analysis/Output/ --no_subspecies
```


#### Taxonomy Searches with and without subspecies:

To understand how the `--no_subspecies` flag can impact analyses, it is important to demonstrate how the taxon list is being parsed. Regardless of the type of names present in this file (species or subspecies), two lists are constructed. One if filled with species (binomial) names, and the other with subspecies (trinomial) names.

Example taxon list file:

```
Leycesteria crocothyrsos
Leycesteria formosa
Linnaea borealis americana
Linnaea borealis borealis
Linnaea borealis longiflora
```

The resulting parsed lists are:

`species = [Leycesteria crocothyrsos, Leycesteria formosa, Linnaea borealis]`

`subspecies = [Linnaea borealis americana, Linnaea borealis borealis, Linnaea borealis longiflora]`

Notice that even though there wasn't a binomial name provided for *Linnaea borealis*, it was automatically generated based on the subspecies labels. This is true regardless of whether the `--no_subspecies` flag is included or not. 

**How does the `--no_subspecies` flag impact searches?**

I will use the above taxon list and the following example records to illustrate:

```
>FJ745393.1 Leycesteria crocothyrsos voucher N. Pyck 1992-1691 maturase K (matK) gene, partial cds; chloroplast
>KC474956.1 Linnaea borealis americana voucher Bennett_06-432_CAN maturase K (matK) gene, partial cds; chloroplast
```

Each sequence record is always parsed to construct a species (binomial) and subspecies (trinomial) name. This would produce the following results:

```
Leycesteria crocothyrsos #species
Leycesteria crocothyrsos voucher #subspecies
Linnaea borealis #species
Linnaea borealis americana #subspecies
```

As you can see above, every subspecies contains a species label. This allows a series of checks to be performed. If the `--no_subspecies` flag is omitted, the following checks are performed:

1. Is the reconstructed species name in the species list? 
    1. If no, the record is ignored.
    2. If yes, the subspecies is examined.
2. Is the reconstructed subspecies name in the subspecies list? 
    1. If no, the species name will be used. 
    2. If yes, the subspecies name will be used instead.

In the example above, `Leycesteria crocothyrsos` is in the species list, but `Leycesteria crocothyrsos voucher` is an obviously incorrect name and is absent from the subspecies list. In this case, the species name `Leycesteria crocothyrsos` will be used for that record. In the other example, `Linnaea borealis` is in the species list, but `Linnaea borealis americana` is also present in the subspecies list, so `Linnaea borealis americana` will be used for that record.

If the `--no_subspecies` flag is included, the following checks are performed:

1. Is the reconstructed species name in the species list? 
    1. If no, the record is ignored.
    2. If yes, the species name is used.

In the example above, `Leycesteria crocothyrsos` and `Linnaea borealis` would be the names used. Essentially, the `--no_subspecies` flag elevates all the subspecies to the species level.

Let's use another example.

Here, the taxon list file contains only species (binomial) names:

```
Draco beccarii
Draco biaro
Draco bimaculatus
Draco blanfordii
Draco boschmai
Draco bourouniensis
Draco cornutus
Draco cristatellus
```

If only binomial names are present, then the resulting species list will be populated and the resulting subspecies list will be empty:

`species = [Draco beccarii, Draco biaro, Draco bimaculatus, Draco blanfordii, Draco boschmai, Draco bourouniensis, Draco cornutus, Draco cristatellus]`

`subspecies = []`

In this example the `--no_subspecies` flag will have no effect on the analysis. That is, regardless of whether the `--no_subspecies` flag is used or not, there aren't any subspecies to reference and the only possible outcome is to find species names.

There are also some special cases depending on combinations of the taxon list and sequence set.

Given the following taxon list:

```
Linnaea borealis americana
Linnaea borealis borealis
Linnaea borealis longiflora
```

And the following record description lines:

```
>KJ593010.1 Linnaea borealis voucher WAB_0132469163 maturase K (matK) gene, partial cds; chloroplast
>KC474956.1 Linnaea borealis americana voucher Bennett_06-432_CAN maturase K (matK) gene, partial cds; chloroplast
>KP297496.1 Linnaea borealis borealis isolate BOP012344 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2, partial sequence
>KP297498.1 Linnaea borealis longiflora isolate BOP022790 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2, partial sequence
```

The following taxa would be detected and included from each record if the `--no_subspecies` flag is omitted:

```
KJ593010.1 -> Linnaea borealis
KC474956.1 -> Linnaea borealis americana
KP297496.1 -> Linnaea borealis borealis
KP297498.1 -> Linnaea borealis longiflora
```

Any record of `Linnaea borealis` missing a valid subspecies label will be lumped in with all other `Linnaea borealis`, whereas those containing valid subspecies labels will be assigned to the correct subspecies. 

The following taxa would be detected and included from each record if the `--no_subspecies` flag is included:

```
KJ593010.1 -> Linnaea borealis
KC474956.1 -> Linnaea borealis
KP297496.1 -> Linnaea borealis
KP297498.1 -> Linnaea borealis
```
This effectively groups all the subspecies under the species name `Linnaea borealis`.

**To summarize:**

+ If the taxon names list contains only species then searches for subspecies labels cannot occur, and therefore the presence or absence of the `--no_subspecies` flag has no effect.
+ If the taxon names list contains a mix of species and subspecies labels, then the `--no_subspecies` flag can substantially change the outcome. 
+ Using the `--no_subspecies` flag reduces all subspecies names to corresponding species names, and is expected to result in less taxa recovered. Depending on your conceptual view of subspecies, you may find this to be an awesome choice, or you may find it to be a terrible choice. 
+ Omitting the `--no_subspecies` flag is expected to produce a greater number of taxa, but only if valid subspecies are present in the starting sequences.
+ There is no downside to having subspecies labels in the taxon list file, because they can effectively be ignored while capturing all relevant species labels.

---------------

### Rename_Merge.py <a name="RM"></a>

The goal of this script is to search through records in a fasta file and replace invalid taxon names with names compatible with the taxon list. The invalid taxon names file (*Unmatched_Taxon_Names.log*) from the `Taxa_Assessment.py` can be used to help create the replacement names file, as it contains all the taxon names that will fail the taxon filter in the `Parse_Loci.py` script. Similar to other input text files, make sure to open the replacement names file in a text editor and ensure the format includes Unix line breaks (line breaks marked by \n, rather than \r\n) and UTF-8 encoding, otherwise extra characters may interfere with parsing the file correctly with **SuperCRUNCH**.

If the optional `-m` flag is used, records that are successfully relabeled are merged with another fasta file. Ideally, this fasta file should be the one containing all the records with valid taxon names (*Matched_Taxa.fasta*). The resulting merged fasta file (*Merged.fasta*) should then be used for the `Parse_Loci.py` script.

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

#### Example Use:

```
python Rename_Merge.py -i bin/Rename/Unmatched_Taxa.fasta -r bin/Rename/taxon_relabeling.txt -o bin/Rename/Output/
```
```
python Rename_Merge.py -i bin/Rename/Unmatched_Taxa.fasta -r bin/Rename/taxon_relabeling.txt -o bin/Rename/Output/ -m bin/Rename/Matched_Taxa.fasta
```


In the replacement names file, the first column should contain the name that needs to be replaced (the invalid name), and the second column should contain the replacement name. Currently, `Rename_Merge.py` only supports species (binomial) name relabeling, so altering subspecies labels is not possible.

Example contents of a replacement names file:

```
Chamaeleo harennae	Trioceros harennae
Chamaeleo hoehneli	Trioceros hoehnelii
Chamaeleo hoehnelii	Trioceros hoehnelii
Chamaeleo jacksonii	Trioceros jacksonii
Chamaeleo johnstoni	Trioceros johnstoni
Chamaeleo melleri	Trioceros melleri
Chamaeleo montium	Trioceros montium
Chamaeleo narraioca	Trioceros narraioca
```

Note that components of a species name are separated with a space, and the columns are separated with a tab character. 

Although the replacement step should help rescue many records, there are some labels in the *Unmatched_Taxon_Names.log* that simply can't be corrected. These include things like:

```
A.alutaceus mitochondrial
A.barbouri mitochondrial
Agama sp.
C.subcristatus tcs1
C.versicolor sox-4
Calotes sp.
Calumma aff.
Calumma cf.
Liolaemus kriegi/ceii
Tsa anolis
Unverified bradypodion
Unverified callisaurus
```

These records have been labeled improperly, or the identity of the organism is uncertain (*sp., cf., aff.*). These are not useful for the analysis, and should rightfully be discarded using the taxonomy filter in `Parse_Loci.py`. Yuck!

In other cases, taxon names may have been updated and now represent synonymies, or may have been accidentally misspelled. Using a organism-specific taxonomy browser can help clarify these situations. These are examples of records that are worth rescuing through relabeling, and using `Rename_Merge.py` to do so will result in higher quality data. 

---------------

### Parse_Loci.py <a name="PL"></a>

To run `Parse_Loci.py`, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms. The goal of this script is to search through the starting sequences and identify records for each gene/marker/locus. Each record that matches a locus must also contain a valid taxon name. 

The taxon names list can contain a mix of species (binomial) and subspecies (trinomial) names. Detailed instructions for the format of this file is provided in the `Taxa_Assessment.py` section. The optional `--no_subspecies` flag can be used, and its effect is also described in great detail in the `Taxa_Assessment.py` section.

The locus file contains the search terms that are used to identify matching records. Detailed information about the format of this file can be found in the *Obtaining Loci and Search Terms* section. 
 
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

#### Example Use:

```
python Parse_Loci.py -i bin/Loci/Merged.fasta -l bin/Loci/locus_search_terms.fasta -t bin/Loci/Taxa_List.txt -o bin/Loci/Output/
```
```
python Parse_Loci.py -i bin/Loci/Merged.fasta -l bin/Loci/locus_search_terms.fasta -t bin/Loci/Taxa_List.txt -o bin/Loci/Output/ --no_subspecies
```


Several output files are created in the directory specified. For each locus included, a fasta file will be written with sequences that pass the locus and taxon filters. If no sequences are found for a locus, a corresponding fasta file will not be produced. A log file summarizing the number of records written per locus is also written to the output directory, and is called *Loci_Record_Counts.log*. An example of the contents of this file is shown below:

```
Locus_Name	Records_Written
BDNF	1246
CMOS	1263
CXCR4	164
EXPH5	650
KIAA1549	0
```
In the example above, the files `BDNF.fasta`, `CMOS.fasta`, `CXCR4.fasta`, and `EXPH5.fasta` will be written to the output directory, but not `KIAA1549.fasta` because no sequences were found. 

Identifying loci in the sequence records through word matching is not expected to be perfect, and there may be non-target sequences in the resulting fasta files. For this reason, the fasta files should be subjected to orthology filtering before attempting to create sequence alignments or performing analyses.

---------------

## **Orthology Filtering** <a name="OF"></a>

![F2](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig2.jpg)

There are two main methods in **SuperCRUNCH** that can be used to filter out paralogous sequences in the locus-specific fasta files. Each relies on using BLAST searches to extract homologous sequence regions, but differ in whether they require user-supplied reference sequences. `Cluster_Blast_Extract.py` works by creating sequence clusters based on similarity using CD-HIT-EST. A BLAST database is constructed from the largest cluster, and all sequences from the fasta file are blasted to this database. In contrast, `Reference_Blast_Extract.py` relies on a previously assembled set of reference sequences to construct the BLAST database, and all sequences from the fasta file are blasted to this database. Both methods offer the ability to specify the BLAST algorithm to use (*blastn, megablast, dc-megablast*), and multiple options for extracting BLAST coordinates. These topics are discussed in greater detail below. `Cluster_Blast_Extract.py` is recommended for 'simple' sequence record sets, in which each record contains sequence data for the same region of a single locus. `Reference_Blast_Extract.py` is recommended for more complex sequence records, in which each record contains multiple loci (long mtDNA fragment, whole organellar genome, etc.) or the records contain non-overlapping fragments of the same locus. The reference sequence set ensures that only the regions of interest are extracted from these records. 

An optional contamination filtering step is also available, called `Contamination_Filter.py`. This step requires a user-supplied set of sequences that represent the source of 'contamination'. For example, amphibian mtDNA sequences can be screened against human mtDNA sequences to make sure they are actually amphibian. Any reference sequences can be used, and the context depends on what the contamination source is expected to be. This step will remove all sequences scoring greater than 95% identity for a minimum of 100 continuous base pairs to the reference ‘contamination’ sequences. 

The resulting filtered fasta files from this step can then be used for the quality filtering and sequence selection stage. 

---------------

### Cluster_Blast_Extract.py <a name="CBE"></a>

`Cluster_Blast_Extract.py` is one of two methods for detecting and removing paralogous sequences from locus-specific fasta files. This method is recommended for 'simple' sequence record sets, in which each record contains sequence data for the same region of a single locus. `Cluster_Blast_Extract.py` works by creating sequence clusters based on similarity using CD-HIT-EST. A BLAST database is constructed from the largest cluster, and all sequences from the fasta file are blasted to this database using the specified BLAST algorithm (*blastn, megablast, dc-megablast*). For each sequence with a significant match, the coordinates of all BLAST hits (excluding self-hits) are merged. This action often results in a single interval, but non-overlapping coordinates can also be produced. Multiple options are provided for handling these cases, with details on this topic provided below. All query sequences with significant hits are extracted based on the resulting BLAST coordinates, and written to a new filtered fasta file. 

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

> **Optional**: The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets. If omitted, no limit is set.

#### Example Use:

```
python Cluster_Blast_Extract.py -i bin/cluster-blast/ -b dc-megablast -m span --max_hits 300
```


Several output folders are created in the directory containing the input fasta files. The directory labels and their contents are described below:

+ **01_Clustering_Results**
    + For each locus, this directory contains the output files from cd-hit-est (ex., `LOCUS1_Out_.clstr`), which delimit the sequence clusters.
+ **02_Parsed_Results**
    + For each locus, this directory contains a set of fasta files which represent the clusters found. These are labeled as `LOCUS1_Out_Cluster_0.fasta`, `LOCUS1_Out_Cluster_1.fasta`, etc.
+ **03_Blast_Results**
    + For each locus, this directory contains the constructed blast database files (extensions .nhr, .nin, .nsq), the blast results for each fasta cluster (ex., `LOCUS1_Out_Cluster_0_blast_results.txt`), and the merged blast results for all clusters (ex., `LOCUS1_blast_results_merged.txt`). 
+ **04_Trimmed_Results**
    + For each locus, this directory contains the filtered fasta file (`LOCUS1_extracted.fasta`) and a corresponding log file (`Log_File_LOCUS1.txt`) that indicates the original sequence length, BLAST coordinates found, and extracted sequence length for each record that passed this filter.
    
#### BLAST Algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use, which can greatly affect the filtering results. The *blastn* algorithm searches with a word size of 11, whereas *megablast* searches include a word size of 28, making blastn more appropriate for interspecies searches and megablast more appropriate for closely related or intraspecific searches. However, discontiguous megablast (*dc-megablast*) is better at producing non-fragmented hits for divergent sequences using similar word sizes as *blastn*, and as such it works well for interspecific and intraspecific searches. If the goal is to produce species level phylogenetic data sets then *dc-megablast* should be used, but if the focus is on population level phylogenetic data sets then *megablast* may be preferable. You can easily compare the effects of the different BLAST algorithms, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

#### BLAST coordinates strategy

The optional `-m ` flag specifies the strategy to use for handling BLAST coordinates. To explain these options, a little background is necessary. When a sequence is blasted, the sections of the query sequence that match the subject sequence (reference) are reported as start and stop coordinates. In general, many hits are produced for the query sequence and the result is many sets of coordinates. Often, these coordinates overlap and can be combined. Take for example the following set of coordinates:

```
[2, 435]
[27, 380]
[30, 500]
```
These can be combined to one set of coordinates: `[2, 500]`. After BLAST searches are finished, the coordinates for each input sequence are merged. Often this produces one set of merged coordinates, like the example above, but sometimes multiple sets of non-overlapping coordinates are produced, like `[2, 70], [100, 500]`. 

How can this happen?

One common reason is the sequence contains a stretch of `N` characters. These will never be matched. Take the following sequence:

```
TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG
```

This sequence will at best produce the following coordinates: `[1, 10], [20, 35]`. The N's are always ignored. 

Another more problematic reason is that duplicate sequences are found, and the coordinates of both duplicates are returned. This can happen in mitochondrial genomes, which sometimes have undergone gene duplications. In cases like this, the genes tend to be separated by quite a long distance. Given a mitochondrial genome, the duplicate genes may return a set of coordinates that look like this: `[300, 860], [4800, 5300]`. The huge gap between these coordinates is a strong signal that duplicate genes are present. 

Given these different scenarios, there are three options meant to handle non-overlapping coordinates, including `span`, `nospan`, and `all`. 

+ When `-m span` is used, if the non-overlapping coordinate sets are less than or equal to 100 base pairs apart, they are merged. In the above 'N' example, this would produce the final coordinates of `[1, 35]`, and the resulting sequence will contain the original stretch of N characters. If the non-overlapping coordinate sets are more than 100 base pairs apart, the coordinates set containing the greatest number of base pairs is selected. In the 'duplication' scenario above, this would produce `[300, 860]`. The `span` option can safely hand the gene duplication scenario, and also allow longer lower-quality sequences to be extracted - as long as they contain a reasonably small stretch of N's. It is the default option used if the `-m ` flag is omitted.

+ When `-m nospan` is used, no attempt is made to merge the non-overlapping coordinates and the coordinate set containing the greatest number of base pairs is selected. In the above 'N' example, this would produce the final coordinates of `[20, 35]` (different from `span`). In the 'duplication' scenario, this would produce `[300, 860]` (same as `span`). The `nospan` option guarantees long stretches of N's will not be present in the extracted sequences, and can penalize these lower-quality sequences by reducing their length. It will also correctly handle the gene duplication scenario. This option should be viewed as a more conservative implementation of `span`.

+ When `-m all` is used, the non-overlapping coordinate sets are used as is. In the 'N' scenario above, this would produce `[1, 10], [20, 35]`, which would simply remove the stretch of N's from the sequence. That is, `TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG` becomes 
`TCATGTTCCACGAAAAATGATGCTG`. For the duplication scenario, this would produce `[300, 860], [4800, 5300]`. In other words, two duplicate genes would be extracted, producing a single sequence that is double the length of normal sequences. For duplications, this is a very poor outcome because it will interfere with sequence alignment. This option can be used to detect duplications in mitogenomes, or paralogous sequences. For example, when I ran this option using many reptile mitogenomes, I was able to find all the gene duplications in mitogenomes for the genus Heteronotia (a parthenogenic gecko). Beyond this use, I would caution against using this option unless you inspect the results very carefully. 

You can easily compare the effects of the options for the `-m ` flag, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

---------------

### Reference_Blast_Extract.py <a name="RBE"></a>

`Reference_Blast_Extract.py` is one of two methods for detecting and removing paralogous sequences from locus-specific fasta files. This method is recommended for more complex sequence records, in which each record contains multiple loci (long mtDNA fragment, whole organellar genome, etc.) or the records contain non-overlapping fragments of the same locus. `Reference_Blast_Extract.py` works by creating a BLAST database from the reference sequences, and all sequences from the fasta file are blasted to this database using the specified BLAST algorithm (*blastn, megablast, dc-megablast*). The reference sequence set ensures that only the target region is extracted from the sequence records. For each sequence with a significant match, the coordinates of all BLAST hits (excluding self-hits) are merged. This action often results in a single interval, but non-overlapping coordinates can also be produced. Multiple options are provided for handling these cases, with details on this topic provided below. All query sequences with significant hits are extracted based on the resulting BLAST coordinates, and written to a new filtered fasta file. 

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

#### Example Use:

```
python Reference_Blast_Extract.py -i bin/Ref-Blast/ -d ND2_references.fasta -e ND2.fasta -b dc-megablast -m span --max_hits 300
```

Several outputs are created in the specified input directory, including:

+ BLAST database files for the reference sequences (extensions .nhr, .nin, .nsq).
+ A BLAST results file, labeled as `LOCUS_blast_output.txt`.
+ The filtered fasta file, labeled as `LOCUS1_extracted.fasta`
+ A corresponding log file (`Log_File_LOCUS1.txt`) that indicates the original sequence length, BLAST coordinates found, and extracted sequence length for each record with significant BLAST hits to the reference sequences.

Similar to `Cluster_Blast_Extract.py`, the same options exist for choosing a BLAST algorithm (`-b `) and the BLAST coordinates strategy (`-m `). For convenience, this information is also posted here. 

#### BLAST Algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use, which can greatly affect the filtering results. The *blastn* algorithm searches with a word size of 11, whereas *megablast* searches include a word size of 28, making blastn more appropriate for interspecies searches and megablast more appropriate for closely related or intraspecific searches. However, discontiguous megablast (*dc-megablast*) is better at producing non-fragmented hits for divergent sequences using similar word sizes as *blastn*, and as such it works well for interspecific and intraspecific searches. If the goal is to produce species level phylogenetic data sets then *dc-megablast* should be used, but if the focus is on population level phylogenetic data sets then *megablast* may be preferable. You can easily compare the effects of the different BLAST algorithms, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

#### BLAST coordinates strategy

The optional `-m ` flag specifies the strategy to use for handling BLAST coordinates. To explain these options, a little background is necessary. When a sequence is blasted, the sections of the query sequence that match the subject sequence (reference) are reported as start and stop coordinates. In general, many hits are produced for the query sequence and the result is many sets of coordinates. Often, these coordinates overlap and can be combined. Take for example the following set of coordinates:

```
[2, 435]
[27, 380]
[30, 500]
```
These can be combined to one set of coordinates: `[2, 500]`
After BLAST, the coordinates for each input sequence are combined. Often this produces one set of merged coordinates, like the above example. Sometimes multiple sets of non-overlapping coordinates are produced, like `[2, 70], [100, 500]`. 

How can this happen?

One common reason is the sequence contains a stretch of `N` characters. These will never be matched. Take the following sequence:

```
TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG
```

This sequence will at best produce the following coordinates: `[1, 10], [20, 35]`. The N's are always ignored. 

Another more problematic reason is that duplicate sequences are found, and the coordinates of both duplicates are returned. This can happen in mitochondrial genomes, which sometimes have undergone gene duplications. In cases like this, the genes tend to be separated by quite a long distance. Given a mitochondrial genome, the duplicate genes may return a set of coordinates that look like this: `[300, 860], [4800, 5300]`. The huge gap between these coordinates is a strong signal that duplicate genes are present. 

Given these different scenarios, there are three options meant to handle non-overlapping coordinates, including `span`, `nospan`, and `all`. 

+ When `-m span` is used, if the non-overlapping coordinate sets are less than or equal to 100 base pairs apart, they are merged. In the above 'N' example, this would produce the final coordinates of `[1, 35]`, and the resulting sequence will contain the original stretch of N characters. If the non-overlapping coordinate sets are more than 100 base pairs apart, the coordinates set containing the greatest number of base pairs is selected. In the 'duplication' scenario above, this would produce `[300, 860]`. The `span` option can safely hand the gene duplication scenario, and also allow longer lower-quality sequences to be extracted - as long as they contain a reasonably small stretch of N's. It is the default option used if the `-m ` flag is omitted.

+ When `-m nospan` is used, no attempt is made to merge the non-overlapping coordinates and the coordinate set containing the greatest number of base pairs is selected. In the above 'N' example, this would produce the final coordinates of `[20, 35]` (different from `span`). In the 'duplication' scenario, this would produce `[300, 860]` (same as `span`). The `nospan` option guarantees long stretches of N's will not be present in the extracted sequences, and can penalize these lower-quality sequences by reducing their length. It will also correctly handle the gene duplication scenario. This option should be viewed as a more conservative implementation of `span`.

+ When `-m all` is used, the non-overlapping coordinate sets are used as is. In the 'N' scenario above, this would produce `[1, 10], [20, 35]`, which would simply remove the stretch of N's from the sequence. That is, `TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG` becomes 
`TCATGTTCCACGAAAAATGATGCTG`. For the duplication scenario, this would produce `[300, 860], [4800, 5300]`. In other words, two duplicate genes would be extracted, producing a single sequence that is double the length of normal sequences. For duplications, this is a very poor outcome because it will interfere with sequence alignment. This option can be used to detect duplications in mitogenomes, or paralogous sequences. For example, when I ran this option using many reptile mitogenomes, I was able to find all the gene duplications in mitogenomes for the genus Heteronotia (a parthenogenic gecko). Beyond this use, I would caution against using this option unless you inspect the results very carefully. 

You can easily compare the effects of the options for the `-m ` flag, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

---------------

### Contamination_Filter.py <a name="CF"></a>

`Contamination_Filter.py` is an optional step for additional filtering, which can help identify and remove 'contaminated' sequences. This step requires a user-supplied set of sequences that represent the source of 'contamination'. For example, amphibian mtDNA sequences can be screened against human mtDNA sequences to make sure they are actually amphibian. Any reference sequences can be used, and the context depends on what the contamination source is expected to be. This step will remove all sequences scoring greater than 95% identity for a minimum of 100 continuous base pairs to the reference ‘contamination’ sequences. 

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

#### Example Use:

```
python Contamination_Filter.py -i bin/contamfilter/ND2/ -d Human_ND2.fasta -e ND2.fasta -b megablast
```

Several outputs are created in the specified input directory, including:

+ BLAST database files for the 'contamination' reference sequences (extensions .nhr, .nin, .nsq).
+ A BLAST results file, labeled as `LOCUS_blast_output.txt`.
+ A filtered fasta file, labeled as `LOCUS1_extracted.fasta`, which contain all sequences that passed the filter.
+ A fasta file of 'contaminated' sequences, labeled as `LOCUS1_extracted_contaminated.fasta`, which contains sequences that failed the filter.
+ A corresponding log file (`Log_File_LOCUS1.txt`), which contains information on the original sequence length, BLAST coordinates found, and extracted sequence length for each record with significant BLAST hits to the reference sequences.

#### BLAST Algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use. For the contamination filter, the goal is to identify and remove sequences with a very high similarity to the references. For this type of search it is best to use *megablast*, which is most appropriate for conducting within-species searches.

---------------

## **Sequence Quality Filtering and Selection** <a name="SQFS"></a>

![F3](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig3.jpg)

The goal of this step is to apply additional sequence filters and select representative sequences using `Filter_Seqs_and_Species.py`. There are multiple options for filtering and selecting sequences, including an option to include or exclude subspecies. Once sequence filtering and selection is completed, the `Make_Acc_Table.py` module can be used to rapidly generate a table of NCBI accession numbers for all taxa and loci.

Although the original goal of `Filter_Seqs_and_Species.py` was to select the best available sequence for a taxon for each locus, it can also be used to retain all sequences passing the filters. That is, it can be used to filter and create population-level data sets. 

---------------

### Filter_Seqs_and_Species.py <a name="FSS"></a>

To construct a phylogenetic supermatrix representative sequences must be selected per taxon per locus, and if multiple sequences exist for a taxon for a locus then an objective strategy must be used for sequence selection. `Filter_Seqs_and_Species.py` includes several strategies for choosing among multiple sequences. The simplest solution is to sort sequences by length and select the longest sequence available (*length* method). `Filter_Seqs_and_Species.py` includes another approach for coding loci in which sequences are first sorted by length, and then undergo translation tests to identify a correct reading frame (*translate* method). Here, the longest sequence is translated in all forward and reverse frames in an attempt to identify a correct reading frame. If a correct reading frame is found, the sequence is selected. If a correct reading frame is not found, the next longest sequence is examined. The process continues until a suitable sequence is found. If no sequences in the set pass translation, the first (longest) sequence is selected (rather than excluding the taxon). The *randomize* feature can be used to shuffle the starting list of sequences randomly, rather than sorting by length. If used in conjunction with the *length* method, the first sequence from the shuffled list is selected (a true random choice). If used in conjunction with the *translate* method, translation tests occur down the shuffled list until a suitable sequence is selected (not necessarily a random choice). 

For all selection options, sequences must meet a minimum base pair threshold set by the user. If no sequences meet the minimum length requirement for a taxon, then the taxon will be excluded during the filtering process.

Similar to previous steps, a list of taxa must be provided. The taxon names list can contain a mix of species (binomial) and subspecies (trinomial) names. Detailed instructions for the format of this file is provided in the `Taxa_Assessment.py` section. The optional `--no_subspecies` flag can be used, and its effect in this step is identical to that described in the `Taxa_Assessment.py` section.

The `Filter_Seqs_and_Species.py` module will create fasta files containing a single representative sequence per taxon for each locus, along with other important output files described below. The filtered fasta files can be used with the `Make_Acc_Table.py` module to generate a table of NCBI accession numbers for all taxa and loci.
 
Although `Filter_Seqs_and_Species.py` was initially designed to filter and select one sequence per species, retaining intraspecific sequence sets may actually be desirable for other projects (e.g., phylogeography). `Filter_Seqs_and_Species.py` includes this option (`--allseqs`). When the `--allseqs` flag is used, rather than selecting the highest quality sequence available for a taxon, all sequences passing the filtration methods are retained. Additional details for this feature are provided below.


#### Basic Usage:

```
python Filter_Seqs_and_Species.py -i <input directory> -f <filter strategy> -l <minimum base pairs> -t <taxon file>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter. The filtering options set are applied to every fasta file in this directory.

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

#### Example Uses:

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --no_subspecies
```
> Above command will select one sequence per taxon per locus based on the longest available sequence. The minimum base pair length is 150. The no_subspecies option is used.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f translate --table standard -l 150 -t bin/Loci/Taxa_List.txt
```
> Above command will select one sequence per taxon per locus based on the longest translatable sequence. The standard code is used for translation for all input fasta files. The minimum base pair length is 150. Subspecies are included.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --no_subspecies --randomize 
```
> Above command will select one sequence per taxon per locus randomly. The minimum base pair length is 150. The no_subspecies option is used.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --allseqs
```
> Above command will select all sequences per taxon per locus that pass the minimum base pair length of 150. Subspecies are included.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f translate --table vertmtdna -l 150 -t bin/Loci/Taxa_List.txt --allseqs
```
> Above command will select all sequences per taxon per locus that pass translation (using vertebrate mitochondrial code) and the minimum base pair length of 150. Subspecies are included.


---------------

### Make_Acc_Table.py <a name="MAT"></a>

Something

#### Basic Usage:

```
python Make_Acc_Table.py -i <input directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the fasta files (single sequence per taxon).

##### `-s <path-to-file>`

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file.





---------------

## **Sequence Alignment** <a name="SA"></a>

![F4](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig4.jpg)

Something


---------------

### Adjust_Direction.py <a name="AD"></a>

Something

#### Basic Usage:

```
python Adjust_Direction.py -i <input directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the unaligned fasta files.

##### `--acc`

> **Optional**: Use the --adjustdirectionaccurately implementation in MAFFT, rather than --adjustdirection. It is slower but more accurate, especially for divergent sequences.


---------------

### Coding_Translation_Tests.py <a name="CTT"></a>

Something

#### Basic Usage:

```
python Coding_Translation_Tests.py -i <input directory> --table <translation table>
```

##### `-i <path-to-directory> `

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter.

##### `--table <choice>`

> **Required**: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*.

##### `--rc`

> **Optional**: In addition to forward frames, use reverse complement for translation tests. Not recommended if direction of sequences has already been adjusted.


---------------

### Align.py <a name="A"></a>

Something

#### Basic Usage:

```
python Align.py -i <input directory> -a <aligner> 
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the unaligned fasta files.

##### `-a <choice>`

> **Required**: Specify whether alignment is by mafft, macse, muscle, or clustalo. If macse, MUST provide flags --mpath and --table. Choices = *mafft, macse, muscle, clustalo, all*.

##### `--mpath <path-to-executable>`

> **Required** for `-a macse`: Full path to a macse jar file.

##### `--table <choice>`

> **Required** for `-a macse`: Specifies translation table. Choices = *standard, vmtdna*.

##### `--mem <integer>`

> **Optional** for `-a macse`: An integer for how much memory to assign to macse (in GB). Default = 1.

##### `--pass_fail`

> **Optional** for `-a macse`: Specifies macse to perform dual alignment. Files in -i directory must follow labeling format: NAME_Passed.fasta, NAME_Failed.fasta.

##### `--accurate`

> **Optional**: Specifies more thorough search settings in mafft, clustalo, or macse (see below).


---------------

### Trim_Alignments.py <a name="TAS"></a>

Something

#### Basic Usage:

```
python Trim_Alignments.py -i <input directory> -f <output format> -a <trimal method>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the input alignment files. File formats can include fasta, nexus, or phylip (see below).

##### `-f <choice>`

> **Required**: Specifies the output file format for trimmed alignments. Choices = *fasta, nexus, phylip*.

##### `-a <choice>`

> **Required**: Specifies the trimal method for trimming alignments. Choices = *gt, noallgaps, both*.

##### `--gt <value>`

> **Optional**: Specifies the gap threshold (gt) value for trimal, the minimum fraction of sequences without a gap. Must be between 0 and 1. Default = 0.05.





---------------

## **File Formatting Tasks** <a name="FFT"></a>


![F5](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig5.jpg)

Something


---------------

### Relabel_Fasta.py <a name="RF"></a>

Something

#### Basic Usage:

```
python Relabel_Fasta.py -i <input directory> -r <relabel option>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the unaligned or aligned fasta files.

##### `-r <choice>`

> **Required**: The strategy for relabeling sequence records. Choices = *species, accession, species_acc*.

##### `-s <path-to-file>`

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file.


---------------

### Fasta_Convert.py <a name="FC"></a>

Something

#### Basic Usage:

```
python Fasta_Convert.py -i <input directory>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the aligned fasta files.


---------------

### Concatenation.py <a name="C"></a>

Something

#### Basic Usage:

```
python Concatenation.py -i <input directory> -r <input format> -s <missing data symbol> -o <output format>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the aligned files.

##### `-f <choice>`

> **Required**: The input file format of the alignments. Choices = *fasta, phylip*.

##### `-s <choice>`

> **Required**: A base pair symbol used to represent missing data sequences. Choices = *dash, N, ?*.

##### `-o <choice>`

> **Required**: The output file format for the final concatenated alignment. Choices = *fasta, phylip*.

---------------

The end.

------
*Last updated: January 2019*
