# Starting Materials

+ [Overview](#GSM)
+ [Obtaining Sequence Data](#OSD)
    + [Using Custom Sequences](#UCS)
    + [Remove_Duplicate_Accessions.py](#RDA)
+ [Obtaining Taxon Names Lists](#OTNL)
    + [Getting Taxa From Fasta Files](#GTFFF)
+ [Obtaining Loci Search Terms](#OLST)
    + [Searching for UCE loci](#SFUL)

---------------

### **Overview** <a name="GSM"></a>

To run **SuperCRUNCH**, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms. Information on how to create these files and recommendations are provided below.

---------------

### Obtaining Sequence Data <a name="OSD"></a>

The starting molecular data for **SuperCRUNCH** consists of a fasta file of nucleotide sequence records. These sequences can be downloaded from the [NCBI nucleotide database (GenBank)](https://www.ncbi.nlm.nih.gov/nucleotide/), generated locally from your own sequencing projects, or can be a combination of both data types. The simplest way to obtain published sequence data is to search for a taxonomic term on the nucleotide database and download all the available records in fasta format. However, for larger taxonomic groups (ex. all birds or frogs) this may not be possible as there will be too much data to download directly. In this case, it is better to split searches using the taxonomy of the group (order, clade, family, etc), download each record set in fasta format, and then combine the resulting fasta files. The [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) is particularly useful for identifying the relevant taxonomic search terms. Advanced users may choose to use additional search terms to filter the results, which could help reduce the resulting file size, but this is not necessary. The downloaded fasta files may be quite large in size (several GB) and it is unlikely you can use a text editor to open and edit them, so combining files requires a different approach. The fasta files will all have the extension '.fasta', and they can be quickly combined using Unix. If the fasta files are all in the same working directory, you can use the following command to accomplish this:

`cat *.fasta > Combined.fasta`

The starting fasta file can be in the interleaved format or sequential format, as demonstrated below: 

**Interleaved format:**
```
>KP820543.1 Callisaurus draconoides voucher MVZ265543 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATTTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTGGGTAATGTAATGGCGAGCGAT
TGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATTCTCAAACAAGAGCAGACAGAGTGCCCCCAGG
AAAATAATTTAATGCTTCCTGAAGACAGCATGGGCATTTTTCAGGATAACAAAAATAATGAACTGTACAA

>KP820544.1 Urosaurus ornatus voucher UWBM7587 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATCTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTAGGTAATGTAATGGCCAGCGAT
TGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATCCTCAAACAAGAACAGACTGAGTGCCCCCAGG
TAAAATTTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTGGGTAATGTAATGGCGAGCGAT
AAACTAATTTAATGCTTCCTGAAGACAGTATGGGTATTTTTCAGGATAACAAAAATAATGAACTGTACAA
```

**Sequential format:**
```
>KP820543.1 Callisaurus draconoides voucher MVZ265543 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATTTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTGGGTAATGTAATGGCGAGCGATTGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATTCTCAAACAAGAGCAGACAGAGTGCCCCCAGGAAAATAATTTAATGCTTCCTGAAGACAGCATGGGCATTTTTCAGGATAACAAAAATAATGAACTGTACAA

>KP820544.1 Urosaurus ornatus voucher UWBM7587 aryl hydrocarbon receptor (anr) gene, partial cds
TAAAATCTCCTTTGAAAGGAACCTTTTTGTGGACACCAGGGATGAATTAGGTAATGTAATGGCCAGCGATTGGCAGGAAAATATTTTGCCAGTGAGGAATAACAGCATCCTCAAACAAGAACAGACTGAGTGCCCCCAGGAAACTAATTTAATGCTTCCTGAAGACAGTATGGGTATTTTTCAGGATAACAAAAATAATGAACTGTACAA
```

The record labels must have a `>` followed by a unique accession number, followed by a sequence description. If curated properly, the description line will contain a taxon name, locus information, and possibly several other identifiers, all separated by whitespace. There should not be any vertical bars (|) separating information contained in the sequence label. **SuperCRUNCH** cannot properly parse this format.

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


### Using Custom Sequences <a name="UCS"></a>

There may be other circumstances in which using a custom set of sequences may be desirable, rather than sequence sets downloaded from NCBI. Any set of sequences can potentially be used in **SuperCRUNCH**. To take advantage of all the available features, the fasta file containing the custom sequences should follow a labeling scheme similar to NCBI sequence records. That is, the general layout for a custom file would involve the following record format:
```
>[Accession] [Genus-label] [Species-label] [optional-Subspecies-label] [Remaining description line (with gene abbreviation and description)]

This is similar to a typical NCBI sequence record:
>KP820543.1 Callisaurus draconoides voucher MVZ265543 aryl hydrocarbon receptor (anr) gene, partial cds
CACAATGAGAAAGCCTTGATAAACCGTGATCGGACTTTGCCACTCGTTGAAGAAATAGATCAGAGGTATT

A custom record could look like this:
>MVZ249370 Hyperolius nitidulus RAG1
CACAATGAGAAAGCCTTGATAAACCGTGATCGGACTTTGCCACTCGTTGAAGAAATAGATCAGAGGTATT
```

Note that all of the above components should be separated by spaces. Because of the way fasta files are read using **SuperCRUNCH**, there cannot be any duplicate accession numbers/identifiers within the same fasta file. If you are using your own custom sequences, make sure to use an appropriate set of identifiers for the accession component. 

The accession component can be any unique identifier using a combination of letters and numbers. In the above example, I used a museum voucher number instead. This would work if the sample is only available for one sequence, but what if the same sample has been sequenced for multiple loci? You can easily create a unique accession number in this case by adding a simple extension to the accession. Let's assume MVZ 249370 has sequence data for 4 loci. The following is an example of how to label the accessions to include all sequences in the starting file that will be used for taxa assessment, parsing loci, etc.

```
>MVZ249370.RAG1 Hyperolius nitidulus recombination activating protein 1 (RAG-1) gene
TTGATAGCTGAAAGAGAAGCCATGAAGTCCAGTGAGCTCATGCTTGAAATCGGCGGAATTCTCAGAAG...

>MVZ249370.TYR Hyperolius nitidulus tyrosinase (TYR) gene
GGGGAAGCCTCAGCTAGAGGAACGTGCCAAGATGTAGTCCTCTCCACTGCTCCTGTAGGTGCTCAATT...

>MVZ249370.POMC Hyperolius nitidulus proopiomelanocortin (POMC) gene
AAAGCATGCAAGATGGACTTATCTGCAGAATCACCTGTATTTCCTGGCAACGGGCACATGCAGCCCCT...

>MVZ249370.FICD Hyperolius nitidulus FIC domain-containing protein (FICD) gene
CACAATGAGAAAGCCTTGATAAACCGTGATCGGACTTTGCCACTCGTTGAAGAAATAGATCAGAGGTA...
```
This labeling will make it clear that all sequences are derived from MVZ 249370, but allow the accessions to be unique. 

The accession component is the only item necessary if you want to run sequences through the orthology filters, but other steps may also require a binomial name (genus label + species label) and/or a description line with loci information to work properly. In the above examples I also provided binomial names and gene abbreviations and descriptions, which would allow these sequences to be incorporated into any step of **SuperCRUNCH**. If custom sequences follow this general format, they can also be merged with NCBI sequences and run together, as highlighted below.
 
For one phylogenetics project, I extracted UCE data from multiple whole genome and also downloaded UCE data from NCBI. Given that there were close to 5,000 loci to parse, filter, and align, I wanted to combine these data gathered from the genome sequences and from NCBI so I could use **SuperCRUNCH** to perform the necessary steps. I labeled the 'custom' sequences from the genomes as follows:

```
>GENOME_Vipera_berus.uce-3517 Vipera berus ultra conserved element uce-3517
TAAGAATGTCGATACCGCTGAACAATCTAAGGGGCTTGGGAAGTGGGGAGCATAAAGTTTAATCTACTTATGCTGGAGTACTTACCTCGT...

>GENOME_Vipera_berus.uce-322 Vipera berus ultra conserved element uce-322
ACAAATATGTGATAAGTAATAAATAAATATCTATTAGTAATCAAAGATCAATAATCAGCCACTGAGCTTACACATTATTCACCATTGTTC...

>GENOME_Vipera_berus.uce-3307 Vipera berus ultra conserved element uce-3307
ACTCACCACTGGCCAATGAGAGTTGAAGTTCATCAACATCTAGTACTTCATTAATTTCCCAACTAAAAGGTATTATGTTTTATGGAAAAT...
```

The main features of these custom records are a unique accession identifier, the genus and species labels, and the appropriate UCE label. Note that for the accession number, I used a period followed by the UCE name to guarantee all the accessions would be unique (similar to the other example above). I did this for all the genomes included. After labeling the custom records in this fashion, I combined them with the typical NCBI records:

```
>KR354177.2 Phrynosoma hernandesi voucher MVZ245875 ultra conserved element uce-1980 genomic sequence
TATCTTCACTGCTACAATTAAAGTAAATTGCAGTTATCTTGTGAGTATCATGTTAGTCTGGTGTTAGAAC
AAACATTTCAATCCTGTTAGGTATAGGACTTCACATTTATTTCAAGTTGCATCACCACATAAATCTG...

>KR354176.2 Phrynosoma modestum voucher MVZ238583 ultra conserved element uce-1980 genomic sequence
TATCCTCACTGCTACAATTAAAGTAAATTGGAGTTATCTTGTGAGTATCATATTAGTTTGGTGTTAGAAC
AAACATTTCAATCCTGTTAGGTATAGGACTTCACATTTATTTCAAGTTGTATCACCACATAAATCTG...

>KR354175.2 Phrynosoma cornutum voucher MVZ238582 ultra conserved element uce-1980 genomic sequence
GTTATCTTGTGAGTATCATGTTAGTCTGGTGTTAGAACAAACATTTCAATCCTGTTAGGTATAGGACTTC
ACATTTATTTCAAGTTGCATCACCACATAAATCTGAGCATGCTAAAATGTATTTCTTTTCCTTATAG...
```

The combined fasta file ran just fine in **SuperCRUNCH**. Note that the custom records were in sequential format, and the NCBI records were in interleaved format - this did not affect the ability of **SuperCRUNCH** to correctly read the fasta file. 

---------------

### Remove_Duplicate_Accessions.py <a name="RDA"></a>

Sometimes combining several fasta files can produce duplicate entries, in other words records with identical accession numbers. This will cause problems with the way **SuperCRUNCH** loads fasta files, as there can be no duplicate accession numbers. You will know this is the case if you try to run a script and it produces the following error:

`ValueError: Duplicate key 'AF443291.2'`

In this instance, there are two records each containing the accession number AF443291.2. In order to use the offending fasta file with **SuperCRUNCH**, you'll need to clean the file to remove any duplicate accession numbers. This script can be used for this purpose.


#### Basic Usage:

```
python Remove_Duplicate_Accessions.py -i <fasta file> -o <output directory>
```

#### Argument Explanations:

##### `-i <full-path-to-file>` 
 
> **Required**:  The full path to a fasta file with GenBank sequence data to filter.

##### `-o <path-to-directory>` 

> **Required**: The full path to an existing directory to write  output fasta file.

#### Example Use:

```
python Remove_Duplicate_Accessions.py -i /bin/starting_material/Iguania_GenBank.fasta -o /bin/starting_material/cleaned/
```
> Above command will remove duplicate sequence records in `Iguania_GenBank.fasta` and the output is written to the specified directory `/bin/starting_material/cleaned/`.

By running this script, all duplicate entries will be removed and a new fasta file will be written to the output directory called *Cleaned.fasta*. This fasta file should be used to run **SuperCRUNCH**.

---------------

### Obtaining Taxon Names Lists <a name="OTNL"></a>

**SuperCRUNCH** requires a list of taxon names that is used to filter sequences. Lists of taxon names can be obtained through general databases, such as the NCBI Taxonomy Browser. In many cases there are specific databases dedicated to major groups, for example the [Reptile Database](http://www.reptile-database.org/), [AmphibiaWeb](https://amphibiaweb.org/), and [Amphibian Species of the World](http://research.amnh.org/vz/herpetology/amphibia/), which usually contain up-to-date taxonomies in a downloadable format. 

The taxon names list required is a simple text file which contains one taxon name per line. The file can contain a mix of species (binomial) and subspecies (trinomial) names, and components of each taxon name should be separated by a space (rather than undescore). Sometimes using Excel or other applications to generate the taxon names text file will include hidden characters not compatible with **SuperCRUNCH**, so make sure to open the file in a text editor and ensure the format includes Unix line breaks (line breaks marked by `\n`, rather than `\r\n`) and UTF-8 encoding. Below are some example contents from suitable taxon name lists.

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

### Getting Taxa From Fasta Files <a name="GTFFF"></a>

In some cases it may be desirable to obtain a list of taxon names directly from a fasta file of sequence records. This option is available using the `Fasta_Get_Taxa.py` module, which is described below. Note that this is inevitably an imperfect process, as messy records can produce spurious names. Although this is a useful shortcut, the resulting list of names should always be inspected and edited before its use in any other steps.

### Fasta_Get_Taxa.py

The goal of this script is to search through all fasta files in a directory to construct 'species' and 'subspecies' label sets directly from the description lines. Two output files are written, which are lists of all unique 'species' and 'subspecies' labels found. 
The names are created by combining the second and third elements on the description line (for binomial names) and by combining the second, third and fourth elements on the description line (for trinomial names). There are several filters in place to try to prevent junk names from being produced, but this is inherently imperfect. Although the filters should work relatively well for constructing binomial names, the trinomial names are a much more difficult problem, and the subspecies list will likely contain some spurious names. The script makes no attempt to address synonomy, and if multiple names are used for the same taxon then they will all be recovered and written to the list.

**The resulting taxon list files should be carefully inspected before using them for any other purpose.** After editing, the lists could be combined to create a taxon list of species and subspecies names. 

For smaller data sets, `Fasta_Get_Taxa.py` can help to generate a taxon list quickly and easily. `Fasta_Get_Taxa.py` was intended to be used for population level data sets, which are unlikely to have a large number of taxa. Smaller size data sets allow for careful inspection for spurious names and synonomies, whereas these task would become arduous at larger taxonomic scales. 


#### Basic Usage:

```
python Fasta_Get_Taxa.py -i <directory with fasta file(s)> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory with fasta file(s) of GenBank sequence data. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write output files.

#### Example Use:

```
python Taxa_Assessment.py -i bin/FastaSet/ -o bin/FastaSet/Output/
```

Two output files are created in the specified output directory, including:

+ `Species_Names.txt`: List of unique binomial names constructed from record descriptions. If records are labeled correctly this should correspond to the genus and species. This file should be inspected.
+ `Subspecies_Names.txt`: List of unique trinomial names constructed from record descriptions. If records actually contain subspecies labels they will be captured in this list, however if the records only contain a binomial name then spurious names may be produced. This file should be VERY carefully inspected.


---------------

### Obtaining Loci and Search Terms <a name="OLST"></a>

**SuperCRUNCH** requires a list of loci and associated search terms to initially identify the content of sequence records. For each locus included in the list, **SuperCRUNCH** will search for associated abbreviated names and longer labels in the sequence record. The choice of loci to include will inherently be group-specific, and surveys of phylogenetic and phylogeographic papers may help to identify an appropriate marker set. There is no limit to the number of loci, and **SuperCRUNCH** can also be used to search for large genomic data sets available on the NCBI nucleotide database, such as those obtained through sequence capture experiments (UCEs, anchored enrichment, etc.). For detailed instructions on searching for UCEs, see the next section.

The format of the locus text file involves three tab-delimited columns. The first column contains the locus name that will be used to label output files. It must not contain any spaces or special characters (including underscores and hyphens), and should be kept short and simple. The second column contains the known abbreviation(s) for the gene or marker. Abbreviations should not include any spaces, but can contain numbers and other characters (like dashes). The second column can contain multiple abbreviations, which should be separated by a semi-colon with no extra spaces between abbreviations. The third column contains a longer label of the gene or marker, such as its full name or description. This third column can also contain multiple search entries, which also should be separated by a semi-colon with no extra spaces between label entries. The abbreviations and labels are not case-specific, as they are converted to uppercase during actual searches along with the description lines of the sequences. Although the locus text file can be created using Excel or other applications, it must be a tab-delimited text file. Similar to the taxon names file, make sure to open the file in a text editor and ensure the format includes Unix line breaks (line breaks marked by `\n`, rather than `\r\n`) and UTF-8 encoding, otherwise extra characters may interfere with parsing the file correctly with **SuperCRUNCH**.

The success of finding loci depends on defining appropriate locus abbreviations and labels. Examples of how searches operate can be found below, which should help guide how to select good search terms. For any locus, it is a good idea to search for the locus on  on GenBank and examine several records to identify the common ways it is labeled.

Here is an example of the formatting for a locus file containing three genes to search for:

```
CMOS	CMOS;C-MOS	oocyte maturation factor
EXPH5	EXPH5	exophilin;exophilin 5;exophilin-5;exophilin protein 5
PTPN	PTPN;PTPN12	protein tyrosine phosphatase;tyrosine phosphatase non-receptor type 12
```

In this example:
+ CMOS contains two abbreviations and one label search term. 
+ EXPH5 contains one abbreviation and four label search terms.
+ PTPN contains two abbreviations and two label search terms. 

The above example is a subset of loci from a locus search terms file I used for squamate reptiles. The complete locus search terms file (`Locus-Search-Terms_Squamate_Markers.txt`) is provided in the example data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data).

**How does the actual locus searching work?**

+ For locus abbreviations, the sequence record label is split by spaces, stripped of punctuation, and converted to uppercase. Each resulting component is checked to see if it is identical to an included locus abbreviation. If so, a match is found. 

+ For locus labels, the sequence record label is converted to uppercase (punctuation is left intact). The line is then checked to see if contains the particular locus label. If so, a match is found. 

+ If a locus abbreviation ***or*** a locus label is matched to the contents of a sequence record, the record will pass the filtering step. 

**Example of locus abbreviation search:**

If the locus file contains:

`CMOS	cmos;c-mos	oocyte maturation factor`

The abbreviations will include:

```
CMOS
C-MOS
```

Notice the search terms are converted to uppercase.

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
Notice the search terms are converted to uppercase.

If the sequence record contains the following label:

`>JX999516.1 Liolaemus pictus voucher LP111; exophilin 5 (EXPH5) gene, partial cds`

It will be converted to the following search line:

`>JX999516.1 LIOLAEMUS PICTUS VOUCHER LP111; EXOPHILIN 5 (EXPH5) GENE, PARTIAL CDS`

Notice punctuation and parentheses are left intact and the line is simply converted to uppercase. The line is then checked to see if any supplied locus label is contained within. In this example, the `EXOPHILIN 5` and `EXOPHILIN` labels are both contained in the line and would produce a match, but `EXOPHILIN-5` and `EXOPHILIN PROTEIN 5` would not. The more specific or complex a label search term is, the less likely it is to produce an exact match. My recommendation is to find the simplest common denominator among records and include that label, along with more complex search labels.

---------------

### Searching for UCE loci <a name="SFUL"></a> 

The strategy for obtaining sets of UCE sequences is a little different from the smaller locus sets. First, you'll want to do a specific GenBank search for the taxonomic term of interest *and* `ultra conserved element` and/or `uce`. This will produce a much more manageable set of sequences to work with, as searching for several thousand loci in a large sequence set will inevitably take a long time. 

To generate a locus search terms file, I retrieved the uce names from the ***uce-5k-probes.fasta*** file located [here](https://github.com/faircloth-lab/uce-probe-sets/tree/master/uce-5k-probe-set). Unfortunately, there does not appear to be a standard naming convention for the UCE loci on GenBank. If the sequences have been properly curated, then the description lines *should* contain the uce name somewhere (uce-10, uce-453, uce-5810, etc). If so, they will be compatible with the the 5k UCE locus search terms file (`Locus-Search-Terms_UCE_5k_set.txt`) I've made available in the data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data).

Here are partial contents from the UCE locus search term file:

```
uce-5805	uce-5805	xxxxxxxxxxxx
uce-5806	uce-5806	xxxxxxxxxxxx
uce-5808	uce-5808	xxxxxxxxxxxx
uce-5810	uce-5810	xxxxxxxxxxxx
```

Notice the third column is junk. Unfortunately, UCE loci have been numbered in a suboptimal way. For example, the label `uce-1` is used instead of `uce-0001`. This causes problems when searching for locus labels, because the term `uce-1` is contained inside of `uce-10`, `uce-104`, `uce-1638`, etc. Because of this, the label search will not work properly, and so we have to rely exclusively on the abbreviation to find the correct sequences.

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

The 5k UCE locus search terms file (`Locus-Search-Terms_UCE_5k_set.txt`) is freely available in the data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data), and it can be used to retrieve UCE data from the 5k set as long as the records have the UCE locus name in the description lines.

---------------