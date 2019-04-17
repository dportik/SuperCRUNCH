# Sequence Quality Filtering and Selection

+ [Overview](#SQFS)
+ [Filter_Seqs_and_Species.py](#FSS)
+ [Make_Acc_Table.py](#MAT)
+ [Infer_Supermatrix_Combinations.py](#ISC)

---------------

## **Overview** <a name="SQFS"></a>

![F3](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig3.jpg)

The goal of this step is to apply additional sequence filters and select representative sequences using `Filter_Seqs_and_Species.py`. There are multiple options for filtering and selecting sequences, including an option to include or exclude subspecies. Once sequence filtering and selection is completed, the `Make_Acc_Table.py` module can be used to rapidly generate a table of NCBI accession numbers for all taxa and loci.

Although the original goal of `Filter_Seqs_and_Species.py` was to select the best available sequence per taxon per locus, it can also be used to retain all sequences passing the filters. That is, it can be used to filter and create population-level data sets, in which all available sequences for each taxon are retained.

If sequences have been filtered for species-level data sets (one sequence per taxon), the total number of possible supermatrix combinations can be calculated using `Infer_Supermatrix_Combinations.py`. This is a fun exercise.

---------------

### Filter_Seqs_and_Species.py <a name="FSS"></a>

To construct a phylogenetic supermatrix representative sequences must be selected per taxon per locus, and if multiple sequences exist for a taxon for a locus then an objective strategy must be used for sequence selection. `Filter_Seqs_and_Species.py` includes several strategies for choosing among multiple sequences. 

+ The simplest solution is to sort sequences by length and select the longest sequence available (`-f length`). 

+ `Filter_Seqs_and_Species.py` includes another approach for coding loci in which sequences are first sorted by length, and then undergo translation tests to identify a correct reading frame (`-f translate`). Here, the longest sequence is translated in all forward and reverse frames in an attempt to identify a correct reading frame. If a correct reading frame is found, the sequence is selected. If a correct reading frame is not found, the next longest sequence is examined. The process continues until a suitable sequence is found. If no sequences in the set pass translation, the first (longest) sequence is selected (rather than excluding the taxon). If the `-f translate` method is used, the translation table must also be specified using the `--table` flag. All NCBI translation table options are available, and can be selected using integers or shortcut terms provided (*standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*). If the `--table` flag is omitted, the default will be to use the Standard translation table. 

+ The `--randomize` feature can be used to shuffle the starting list of sequences randomly, rather than sorting by length. If used in conjunction with the `-f length` method, the first sequence from the shuffled list is selected (a true random choice). If used in conjunction with the `-f translate` method, translation tests occur down the shuffled list until a suitable sequence is selected (not necessarily a random choice). 

For all selection options, sequences must meet a minimum base pair threshold set by the user. If no sequences meet the minimum length requirement for a taxon, then the taxon will be eliminated during the filtering process. To avoid eliminating short sequences a small integer can be used, although the retention of very short sequences can interfere with sequence alignment.

Similar to previous steps, a list of taxa must be provided. The taxon names list can contain a mix of species (binomial) and subspecies (trinomial) names. Detailed instructions for the format of this file is provided in the `Taxa_Assessment.py` section. The optional `--no_subspecies` flag can be used, and its effect in this step is identical to that described in the `Taxa_Assessment.py` section.

Using one of the sequence selection strategies above, the `Filter_Seqs_and_Species.py` module will create fasta files containing a single representative sequence per taxon for each locus, along with other important output files described below. The filtered fasta files can be used with the `Make_Acc_Table.py` module to generate a table of NCBI accession numbers for all taxa and loci.
 
Although `Filter_Seqs_and_Species.py` was initially designed to filter and select one sequence per species, retaining intraspecific sequence sets may actually be desirable for other projects (e.g., phylogeography). `Filter_Seqs_and_Species.py` includes the option to create these data sets using the `--allseqs` flag. When the `--allseqs` flag is used, rather than selecting the highest quality sequence available for a taxon, all sequences passing the filtration methods are retained. 


#### Basic Usage:

```
python Filter_Seqs_and_Species.py -i <input directory> -f <filter strategy> -l <minimum base pairs> -t <taxon file>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter. The filtering options set are applied to every fasta file in this directory. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `-f <choice>`

> **Required**: Strategy for filtering sequence data. Choices = *translate, length*.

##### `-l <integer>`

> **Required**: An integer for the minimum number of base pairs required to keep a sequence (ex. 150).

##### `-t <path-to-file>`

> **Required**: The full path to a text file containing all taxon names to cross-reference in the fasta file.

##### `--table <choice>`

> **Required for** `-f translate`: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*.

##### `--randomize`

> **Optional**: For taxa with multiple sequences, shuffle order randomly. Overrides sorting by length for all methods specified using `-f`.

##### `--allseqs`

> **Optional**: For taxa with multiple sequences, select all sequences passing the filters instead of a single representative sequence.

##### `--no_subspecies`

> **Optional**: Ignore subspecies labels in both the taxon names file and the fasta file.

#### Example Uses:

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --no_subspecies
```
> Above command will select one sequence per taxon per locus based on the longest available sequence. The minimum base pair length is 150. The subspecies component of taxon names is ignored. Action is performed for all fasta files located in the `Filter/` directory.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f translate --table standard -l 150 -t bin/Loci/Taxa_List.txt
```
> Above command will select one sequence per taxon per locus based on the longest translatable sequence. The standard code is used for translation for all input fasta files. The minimum base pair length is 150. Subspecies are included. Action is performed for all fasta files located in the `Filter/` directory.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --no_subspecies --randomize 
```
> Above command will select one sequence per taxon per locus randomly. The minimum base pair length is 150. The subspecies component of taxon names is ignored. Action is performed for all fasta files located in the `Filter/` directory.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f length -l 150 -t bin/Loci/Taxa_List.txt --allseqs
```
> Above command will select **ALL** sequences per taxon per locus that pass the minimum base pair length of 150. Subspecies are included. Action is performed for all fasta files located in the `Filter/` directory.

```
python Filter_Seqs_and_Species.py -i /bin/Filter/ -f translate --table vertmtdna -l 150 -t bin/Loci/Taxa_List.txt --allseqs
```
> Above command will select **ALL** sequences per taxon per locus that pass translation (using vertebrate mitochondrial code) and the minimum base pair length of 150. Subspecies are included. Action is performed for all fasta files located in the `Filter/` directory.

Several outputs are created in the specified input directory (one for every input fasta file):

+ `[fasta name]_single_taxon.fasta`: An output fasta file which contains a single filtered sequence per taxon. Produced if the `--allseqs` flag **is not used**.

**OR**

+ `[fasta name]_all_seqs.fasta`: An output fasta file which contains all available sequences per taxon that passed the relevant filters. Produced if the `--allseqs` flag **is used**.

+ `[fasta name]_species_log.txt`: A summary file containing information for each taxon entry. An example of the contents is shown below. In addition to reporting the sequence accession number and length, the number of alternative sequences available is shown. If the `-f translate` method was used, the results of whether the translation test was passed is provided (Y or N). If the `-f length` method was used, this column will contain NA instead.
```
Taxon	Accession	SeqLength	PassedTranslation	SeqsAvailable
Agama lionotus	GQ242168.1	853	Y	1
Amblyrhynchus cristatus	NC_028031.1	1038	Y	21
Amphibolurus muricatus	HQ684202.1	1032	Y	66
Amphibolurus norrisi	AY133001.1	1032	Y	18
Anisolepis longicauda	AF528736.1	1033	Y	1
Anolis acutus	AF055926.2	1029	Y	1
Anolis aeneus	AF055950.1	1033	Y	2
...
```

+ `[fasta name]_accession_list_by_species.txt`: A tab-delimited file in which each line starts with a taxon name and is followed by all accession numbers of sequences passing the length filter from the fasta file. These are the accession numbers for the selected sequence and all alternative sequences indicated in the file above. An example of the contents is shown below:
```
Agama lionotus	GQ242168.1	
Amblyrhynchus cristatus	NC_028031.1	KT277937.1	KR350765.1	KR350766.1	KR350762.1	KR350754.1	KR350759.1	KR350763.1	KR350764.1	KR350758.1	KR350755.1	KR350761.1	KR350753.1	KR350757.1	KR350752.1	KR350760.1	KR350756.1	KR350743.1	KR350746.1	KR350744.1	KR350745.1	
Amphibolurus muricatus	HQ684202.1	HQ684207.1	HQ684201.1	AF128468.1	HQ684203.1	HQ684199.1	HQ684206.1	HQ684204.1	HQ684205.1	KF871666.1	KF871689.1	KF871681.1	KF871658.1	KF871665.1	KF871679.1	KF871692.1	KF871706.1	KF871655.1	KF871688.1	KF871659.1	KF871656.1	KF871653.1	KF871670.1	KF871703.1	KF871699.1	KF871691.1	KF871685.1	KF871701.1	KF871669.1	KF871675.1	KF871705.1	KF871704.1	KF871696.1	KF871700.1	KF871671.1	KF871676.1	KF871702.1	KF871672.1	KF871660.1	KF871668.1	KF871697.1	KF871683.1	KF871667.1	KF871677.1	KF871687.1	KF871698.1	KF871686.1	KF871678.1	KF871664.1	KF871694.1	KF871651.1	KF871684.1	KF871662.1	HQ684200.1	KF871654.1	KF871661.1	KF871680.1	KF871693.1	KF871663.1	KF871673.1	KF871690.1	KF871652.1	KF871650.1	KF871657.1	KF871682.1	KF871695.1	
Amphibolurus norrisi	AY133001.1	HQ684197.1	HQ684208.1	HQ684198.1	HQ684211.1	HQ684210.1	HQ684195.1	HQ684194.1	HQ684188.1	HQ684196.1	HQ684192.1	HQ684191.1	HQ684190.1	HQ684189.1	HQ684209.1	HQ684193.1	KF871674.1	HQ684187.1	
Anisolepis longicauda	AF528736.1	
Anolis acutus	AF055926.2	
Anolis aeneus	AF055950.1	AF317066.1	
...
```

+ `[fasta name]_accession_list_for_Batch_Entrez.txt`: a Batch Entrez style file which contains all the accession numbers for sequences that passed the length filter. This list is a combination of all the accession numbers in the above file. It can be used to download all records using the Batch Entrez portal on NCBI. An example of the contents is shown below:
```
GQ242168.1
NC_028031.1
KT277937.1
KR350765.1
KR350766.1
KR350762.1
KR350754.1
KR350759.1
KR350763.1
KR350764.1
KR350758.1
...
```

These output files provide explicit information regarding which sequence has been selected for all taxa and loci. In addition, they provide accession numbers that can be used to re-create the fasta files used.

---------------

### Make_Acc_Table.py <a name="MAT"></a>

The goal of `Make_Acc_Table.py` is to create a table of accession numbers for each taxon and locus from a directory of fasta files, where each fasta file represents a different locus. These files must have been previously screened such that each taxon is only represented by a single sequence per locus. 

A complete set of taxa is inferred from all the fasta files. By default the names constructed are binomial (genus and species). If subspecies names are also desired, the optional `-s ` flag must be used and a text file of subspecies names must be supplied. The same taxon list from previous steps can be used here, but only the subspecies names are extracted from this list. If a record matches an included subspecies label then the subspecies name will be used, otherwise the binomial name will be used. 

For every taxon, the accession numbers are obtained from each locus. If the taxon is not present in a fasta file (locus), it will appear in the table as a dash. The names of the columns will match the fasta file names, and the rows are composed of taxon names in alphabetical order.

**Note**: This script can process unaligned or aligned fasta files, but they must have the original description lines (not relabeled), and taxa must only contain a single sequence per locus.

#### Basic Usage:

```
python Make_Acc_Table.py -i <input directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the fasta files (single sequence per taxon). Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `-s <path-to-file>`

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file. This can be a taxon names file used in previous steps, or a smaller version only containing subspecies names.

#### Example Uses:

```
python Make_Acc_Table.py -i /bin/filtered_fasta/ -s /bin/taxa/taxon_list.txt
```
> Above command will construct a table of accession numbers based on the fasta files present in the `filtered_fasta/` directory. Subspecies labels from `taxon_list.txt` will be used to find subspecies names in the fasta files.

The accession table will be written as a tab-delimited text file to the input directory as `GenBank_Accession_Table.txt`. Here is an example of the contents from a truncated file:

```
Taxon	12S	16S	CO1	CYTB	ND1
Acanthocercus adramitanus	-	KU097508.1	-	-	-
Acanthocercus annectens	-	MG700133.1	MG699914.1	-	-
Acanthocercus atricollis	-	JX668132.1	-	-	-
Acanthocercus cyanogaster	-	JX668138.1	-	-	-
Acanthocercus yemensis	-	JX668140.1	-	-	-
Acanthosaura armata	NC_014175.1	NC_014175.1	NC_014175.1	NC_014175.1	NC_014175.1
Acanthosaura capra	-	-	-	AY572880.1	-
Acanthosaura crucigera	AB031963.1	MG935713.1	MG935416.1	AY572889.1	-
Acanthosaura lepidogaster	KR092427.1	KR092427.1	KR092427.1	KR092427.1	KR092427.1
Agama aculeata	-	JX668143.1	-	AF355563.1	-
...
```

This file can be opened and manipulated using other applications such as Excel. By default, the columns are sorted in alphabetical order according to the fasta names.


---------------

### Infer_Supermatrix_Combinations.py <a name="ISC"></a>

The goal of `Infer_Supermatrix_Combinations.py` is to calculate how many supermatrix combinations are available, given the number of filtered sequences available for each taxon for each locus. If all taxa have only one sequence available, the answer is one, but if taxa have multiple sequences available, this number will be extremely large. This module relies on the `[locus]_species_log.txt` files produced from the `Filter_Seqs_and_Species.py` module to calculate the number of sequences available per taxon. The log files for all loci should be present in the input directory for the calculation to be accurate.

No output files are created, rather the information is logged to the screen. This includes the total number of sequences available (for all taxa across all loci), the total number of taxa, and the total number of possible supermatrix combinations.

#### Basic Usage:

```
python Infer_Supermatrix_Combinations.py -i <input directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains all the `[locus]_species_log.txt` files.

#### Example Uses:

```
python Infer_Supermatrix_Combinations.py -i /bin/filtered_files/
```
> Above command will calculate the total number of supermatrix combinations based on the  `[locus]_species_log.txt` files present in the `filtered_files/` directory.

Example output on screen:

```
Found 4 loci to examine.


	Parsing information in ITS_extracted_species_log.txt
	Parsing information in MATK_extracted_species_log.txt
	Parsing information in RBCL_extracted_species_log.txt
	Parsing information in TRNL-TRNF_extracted_species_log.txt


Found 4,098 total sequences for 651 taxa.


Number of possible supermatrix combinations (unwieldy integer) = 23,585,393,330,101,509,977,748,274,541,526,924,398,183,623,226,243,709,268,339,387,872,606,971,131,321,031,032,654,793,314,825,517,791,632,690,451,488,529,178,270,403,709,016,922,288,924,652,228,738,591,235,421,706,300,132,306,309,071,552,868,337,445,042,676,361,861,581,482,286,906,518,734,986,582,600,423,724,641,187,945,561,979,141,257,183,422,301,145,122,812,952,763,816,817,181,582,943,498,141,376,096,275,064,842,431,024,332,800,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000.


Number of possible supermatrix combinations (scientific notation) = 2.36E+391, or 2.36*10^391.


```

That's a lot of possible supermatrices! Good thing **SuperCRUNCH** has an objective, repeatable method for selecting sequences.

---------------
