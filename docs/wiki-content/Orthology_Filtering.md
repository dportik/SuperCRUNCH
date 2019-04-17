# Orthology Filtering

+ [Overview](#OF)
+ [Cluster_Blast_Extract.py](#CBE)
+ [Reference_Blast_Extract.py](#RBE)
+ [Contamination_Filter.py](#CF)

---------------

## **Overview** <a name="OF"></a>

![F2](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig2.jpg)

There are two main methods in **SuperCRUNCH** that can be used to filter out non-orthologous sequences in the locus-specific fasta files. Each relies on using BLAST searches to identify and extract homologous sequence regions, but they differ in whether they require user-supplied reference sequences. `Cluster_Blast_Extract.py` works by creating sequence clusters based on similarity using CD-HIT-EST. A BLAST database is automatically constructed from the largest cluster, and all sequences from the fasta file are blasted to this database. In contrast, `Reference_Blast_Extract.py` relies on a previously assembled set of reference sequences to construct the BLAST database, and all sequences from the fasta file are blasted to this database. Both methods offer the ability to specify the BLAST algorithm to use (*blastn, megablast, dc-megablast*), and multiple options for extracting sequence targets based on resulting BLAST coordinates. These topics are discussed in greater detail below. `Cluster_Blast_Extract.py` is recommended for 'simple' sequence record sets, in which each record contains sequence data for the same region of a single locus. `Reference_Blast_Extract.py` is recommended for 'complex' sequence records, in which records contain multiple loci (long mtDNA fragments, whole organellar genome, etc.), non-overlapping fragments of the same locus (from different primer sets), or considerable variation in the length of the target region across records. The reference sequence set ensures that only the target regions are extracted from these records. Examples of reference sequence sets are available in the data folder [here](https://github.com/dportik/SuperCRUNCH/tree/master/data).

An optional contamination filtering step is also available, called `Contamination_Filter.py`. This step requires a user-supplied set of sequences that represent the source of 'contamination'. For example, amphibian mtDNA sequences can be screened against human mtDNA sequences to make sure they are actually amphibian. Any set of reference sequences can be used, and the context depends on what the contamination source is expected to be. This step will remove all sequences scoring greater than 95% identity for a minimum of 100 continuous base pairs to the reference ‘contamination’ sequences. 

The fasta files resulting from this step can be used for the quality filtering and sequence selection stage. 

---------------

### Cluster_Blast_Extract.py <a name="CBE"></a>

`Cluster_Blast_Extract.py` is one of two methods for detecting and removing non-orthologous sequences from locus-specific fasta files. This method is recommended for simple sequence record sets, in which each record contains sequence data for the same region of a single locus. `Cluster_Blast_Extract.py` works by creating sequence clusters based on similarity using CD-HIT-EST. A BLAST database is constructed from the largest cluster, and all sequences from the fasta file are blasted to this database using the specified BLAST algorithm (*blastn, megablast,* or *dc-megablast*). For each sequence with a significant match, the coordinates of all BLAST hits (excluding self-hits) are merged. This action often results in a single interval, but non-overlapping coordinates can also be produced. Multiple options are provided for handling these cases, with details on this topic provided below. All query sequences with significant hits are extracted based on the resulting BLAST coordinates, and written to a new filtered fasta file. If any sequences failed this step, they are written to a separate output fasta file.

#### Basic Usage:

```
python Cluster_Blast_Extract.py -i <fasta file directory> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the parsed locus-specific fasta files. Fasta files in the directory must have the extension '.fasta' to be read.

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. **Recommended**: *dc-megablast*.

##### `-m <choice>`

> **Optional**: The strategy for dealing with multiple non-overlapping blast coordinates. Choices = *span, nospan, all*. Default = *span*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets. If omitted, no limit is set.

#### Example Use:

```
python Cluster_Blast_Extract.py -i bin/cluster-blast/ -b dc-megablast -m span --max_hits 300
```
> Above command will perform automated clustering, BLASTing using *dc-megablast* (with hit limit imposed), and the *span* strategy for BLAST coordinates for each unaligned fasta file present in the `cluster-blast/` input directory.


Several output folders are created in the directory containing the input fasta files. The directory labels and their contents are described below:

+ **01_Clustering_Results/**
    + For each locus, this directory contains the output files from cd-hit-est (ex., `LOCUS1_Out.clstr`), which delimit the sequence clusters.
+ **02_Parsed_Results/**
    + For each locus, this directory contains a set of fasta files which represent the clusters found. These are labeled as `LOCUS1_Out_Cluster_0.fasta`, `LOCUS1_Out_Cluster_1.fasta`, etc.
+ **03_Blast_Results/**
    + For each locus, this directory contains the constructed blast database files (extensions .nhr, .nin, .nsq), the blast results for each fasta cluster (ex., `LOCUS1_Out_Cluster_0_blast_results.txt`), and the merged blast results for all clusters (ex., `LOCUS1_blast_results_merged.txt`). 
+ **04_Trimmed_Results/**
    + For each locus, this directory contains the filtered fasta file (`LOCUS1_extracted.fasta`) and a corresponding log file (`Log_File_LOCUS1.txt`) that indicates the original sequence length, BLAST coordinates found, and extracted sequence length for each record that passed this filter. Additionally, if any sequences fail the orthology filter an additional fasta file is created (`Log_BadSeqs_LOCUS1.fasta`). These records can be inspected to identify the potential reasons for not passing the filter.
    
#### BLAST algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use, which can greatly affect the filtering results. The *blastn* algorithm searches with a word size of 11, whereas *megablast* searches include a word size of 28, making *blastn* more appropriate for interspecies searches and *megablast* more appropriate for closely related or intraspecific searches. However, *discontiguous megablast* (*dc-megablast*) is better at producing non-fragmented hits for divergent sequences using similar word sizes as *blastn*, and as such it works well for interspecific and intraspecific searches. If the goal is to produce species level phylogenetic data sets then *dc-megablast* should be used, but if the focus is on population level phylogenetic data sets then *megablast* may be preferable. You can easily compare the effects of the different BLAST algorithms, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

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

+ When `-m span` is used, if the non-overlapping coordinate sets are less than or equal to 100 base pairs apart, they are merged. In the above 'N' example, this would produce the final coordinates of `[1, 35]`, and the resulting sequence will contain the original stretch of N characters. If the non-overlapping coordinate sets are more than 100 base pairs apart, the coordinates set containing the greatest number of base pairs is selected. In the 'duplication' scenario above, this would produce `[300, 860]`. The `span` option can safely handle the gene duplication scenario, and also allow longer lower-quality sequences to be extracted - as long as they contain a reasonably small stretch of N's. It is the default option used if the `-m ` flag is omitted.

+ When `-m nospan` is used, no attempt is made to merge the non-overlapping coordinates and the coordinate set containing the greatest number of base pairs is selected. In the above 'N' example, this would produce the final coordinates of `[20, 35]` (different from `span`). In the 'duplication' scenario, this would produce `[300, 860]` (same as `span`). The `nospan` option guarantees long stretches of N's will not be present in the extracted sequences, and can penalize these lower-quality sequences by reducing their length. It will also correctly handle the gene duplication scenario. This option should be viewed as a more conservative implementation of `span`.

+ When `-m all` is used, the non-overlapping coordinate sets are used as is. In the 'N' scenario above, this would produce `[1, 10], [20, 35]`, which would simply remove the stretch of N's from the sequence. That is, `TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG` becomes 
`TCATGTTCCACGAAAAATGATGCTG`. Although this seems like a desirable outcome, the same strategy will severely backfire for duplications. For the duplication scenario, this would produce `[300, 860], [4800, 5300]`. In other words, two duplicate genes would be extracted, producing a single sequence that is double the length of normal sequences. For duplications, this is a very poor outcome because it will interfere with sequence alignment. This option can be used to detect duplications in mitogenomes, or paralogous sequences. For example, when I ran this option using many reptile mitogenomes, I was able to find all the gene duplications in mitogenomes for the genus *Heteronotia* (a parthenogenic gecko). Beyond this use, I would caution against using this option unless you inspect the results very carefully. 

You can easily compare the effects of the options for the `-m ` flag, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

---------------

### Reference_Blast_Extract.py <a name="RBE"></a>

`Reference_Blast_Extract.py` is one of two methods for detecting and removing non-orthologous sequences from locus-specific fasta files. This method is recommended for more complex sequence records, in which some records contain multiple loci (long mtDNA fragment, whole organellar genome, etc.) or the records contain non-overlapping fragments of the same locus. `Reference_Blast_Extract.py` works by creating a BLAST database from the reference sequences, and all sequences from the fasta file are blasted to this database using the specified BLAST algorithm (*blastn, megablast,* or *dc-megablast*). The reference sequence set ensures that only the target region is extracted from the sequence records. For each sequence with a significant match, the coordinates of all BLAST hits (excluding self-hits) are merged. This action often results in a single interval, but non-overlapping coordinates can also be produced. Multiple options are provided for handling these cases, with details on this topic provided below. All query sequences with significant hits are extracted based on the resulting BLAST coordinates, and written to a new filtered fasta file. If any sequences failed this step, they are written to a separate output fasta file.

#### Basic Usage:

```
python Reference_Blast_Extract.py -i <input directory> -d <reference fasta name> -e <empirical fasta name> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the reference fasta file and the empirical fasta file.

##### `-d <filename>`

> **Required**: The name of the reference fasta file that will be used to create the blast database. Requires file name only. **DO NOT PROVIDE A FULL FILE PATH**. This file should be located in the input directory (-i), and it will identified using only the file name.

##### `-e <filename>`

> **Required**: The name of the empirical fasta file to blast to the database to prune sequences. Requires file name only. **DO NOT PROVIDE A FULL FILE PATH**. This file should be located in the input directory (-i), and it will identified using only the file name.

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. **Recommended**: *dc-megablast*.

##### `-m <choice>`

> **Optional**: The strategy for dealing with multiple non-overlapping blast coordinates. Choices = *span, nospan, all*. Default = *span*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence. May want to set < 300 for large sequence sets.

#### Example Use:

```
python Reference_Blast_Extract.py -i bin/Ref-Blast/ -d ND2_references.fasta -e ND2.fasta -b dc-megablast -m span --max_hits 300
```
> Above command will form a BLAST database from `ND2_references.fasta` and BLAST sequences from `ND2.fasta` to the database using the *dc-megablast* algorithm, *span* strategy for BLAST coordinate merging, and a hit limit of 300. Both files are in the `bin/Ref-Blast/` directory, where outputs are also written.

Several outputs are created in the specified input directory, including:

+ BLAST database files for the reference sequences (extensions .nhr, .nin, .nsq).
+ A BLAST results file, labeled as `[LOCUS1]_blast_output.txt`.
+ The filtered fasta file, labeled as `[LOCUS1]_extracted.fasta`
+ A corresponding log file (`Log_File_[LOCUS1].txt`) that indicates the original sequence length, BLAST coordinates found, and extracted sequence length for each record with significant BLAST hits to the reference sequences.
+ Additionally, if any sequences fail the orthology filter an additional fasta file is created (`Log_BadSeqs_[LOCUS1].fasta`). These records can be inspected to identify the potential reasons for not passing the filter.

Similar to `Cluster_Blast_Extract.py`, the same options exist for choosing a BLAST algorithm (`-b `) and the BLAST coordinates strategy (`-m `). For convenience, this information is also posted here. 

#### BLAST algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use, which can greatly affect the filtering results. The *blastn* algorithm searches with a word size of 11, whereas *megablast* searches include a word size of 28, making *blastn* more appropriate for interspecies searches and *megablast* more appropriate for closely related or intraspecific searches. However, *discontiguous megablast* (*dc-megablast*) is better at producing non-fragmented hits for divergent sequences using similar word sizes as *blastn*, and as such it works well for interspecific and intraspecific searches. If the goal is to produce species level phylogenetic data sets then *dc-megablast* should be used, but if the focus is on population level phylogenetic data sets then *megablast* may be preferable. You can easily compare the effects of the different BLAST algorithms, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

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

+ When `-m span` is used, if the non-overlapping coordinate sets are less than or equal to 100 base pairs apart, they are merged. In the above 'N' example, this would produce the final coordinates of `[1, 35]`, and the resulting sequence will contain the original stretch of N characters. If the non-overlapping coordinate sets are more than 100 base pairs apart, the coordinates set containing the greatest number of base pairs is selected. In the 'duplication' scenario above, this would produce `[300, 860]`. The `span` option can safely handle the gene duplication scenario, and also allow longer lower-quality sequences to be extracted - as long as they contain a reasonably small stretch of N's. It is the default option used if the `-m ` flag is omitted.

+ When `-m nospan` is used, no attempt is made to merge the non-overlapping coordinates and the coordinate set containing the greatest number of base pairs is selected. In the above 'N' example, this would produce the final coordinates of `[20, 35]` (different from `span`). In the 'duplication' scenario, this would produce `[300, 860]` (same as `span`). The `nospan` option guarantees long stretches of N's will not be present in the extracted sequences, and can penalize these lower-quality sequences by reducing their length. It will also correctly handle the gene duplication scenario. This option should be viewed as a more conservative implementation of `span`.

+ When `-m all` is used, the non-overlapping coordinate sets are used as is. In the 'N' scenario above, this would produce `[1, 10], [20, 35]`, which would simply remove the stretch of N's from the sequence. That is, `TCATGTTCCANNNNNNNNNNCGAAAAATGATGCTG` becomes 
`TCATGTTCCACGAAAAATGATGCTG`. Although this seems like a desirable outcome, the same strategy will severely backfire for duplications. For the duplication scenario, this would produce `[300, 860], [4800, 5300]`. In other words, two duplicate genes would be extracted, producing a single sequence that is double the length of normal sequences. For duplications, this is a very poor outcome because it will interfere with sequence alignment. This option can be used to detect duplications in mitogenomes, or paralogous sequences. For example, when I ran this option using many reptile mitogenomes, I was able to find all the gene duplications in mitogenomes for the genus *Heteronotia* (a parthenogenic gecko). Beyond this use, I would caution against using this option unless you inspect the results very carefully. 

You can easily compare the effects of the options for the `-m ` flag, as the coordinates used to extract sequences are readily available in the log files produced in the final `/04_Trimmed_Results` directory.

---------------

### Contamination_Filter.py <a name="CF"></a>

`Contamination_Filter.py` is an optional step for additional filtering, which can help identify and remove 'contaminated' sequences. This step requires a user-supplied set of sequences that represent the source of 'contamination'. For example, amphibian mtDNA sequences can be screened against human mtDNA sequences to make sure they are actually amphibian. Any reference sequences can be used, and the context depends on what the contamination source is expected to be. This step will remove all sequences scoring greater than 95% identity for a minimum of 100 continuous base pairs to the reference ‘contamination’ sequences. 

In general, it is easier to create one database of the 'contamination source' that contains all the genes you will want to BLAST against. For the Iguania data set, I used a contamination reference composed of human mtDNA, with one sequence for 12S, 16S, CO1, CYTB, ND1, ND2, and ND4. This file `Contamination_Seqs_human_mtDNA.fasta` is available in the [data folder](https://github.com/dportik/SuperCRUNCH/tree/master/data). 

#### Basic Usage:

```
python Contamination_Filter.py -i <input directory> -d <contamination fasta name> -e <empirical fasta name> -b <blast algorithm> -m <blast coordinate strategy>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the reference fasta file and the empirical fasta file.

##### `-d <filename>`

> **Required**: The name of the contamination reference fasta file that will be used to create the blast database. Requires file name only. **DO NOT PROVIDE A FULL FILE PATH**. This file should be located in the input directory (-i), and it will identified using only the file name.

##### `-e <filename>`

> **Required**: The name of the empirical fasta file to blast to the database to find bad sequences. Requires file name only. **DO NOT PROVIDE A FULL FILE PATH**. This file should be located in the input directory (-i), and it will identified using only the file name.

##### `-b <choice>`

> **Required**: The blast algorithm to use. Choices = *blastn, blastn-short, dc-megablast, megablast*. **Recommended**: *megablast*.

##### `--max_hits <integer>`

> **Optional**: The maximum number of blast matches allowed per input sequence.

#### Example Use:

```
python Contamination_Filter.py -i bin/contamfilter/ND2/ -d Human_ND2.fasta -e ND2.fasta -b megablast
```
> Above command will form a BLAST database from `Human_ND2.fasta` and BLAST sequences from `ND2.fasta` to the database using the *megablast* algorithm. Both files are in the `bin/Ref-Blast/` directory, where output files are also written.

Several outputs are created in the specified input directory, including:

+ BLAST database files for the 'contamination' reference sequences (extensions .nhr, .nin, .nsq).
+ A BLAST results file, labeled as `[LOCUS1]_blast_output.txt`.
+ A filtered fasta file, labeled as `[LOCUS1]_extracted.fasta`, which contain all sequences that passed the filter.
+ A fasta file of 'contaminated' sequences, labeled as `[LOCUS1]_extracted_contaminated.fasta`, which contains sequences that failed the filter.
+ A corresponding log file (`Log_File_[LOCUS1].txt`), which contains information on the original sequence length, BLAST coordinates found, and extracted sequence length for each record with significant BLAST hits to the reference sequences.

#### BLAST algorithm choice

The required `-b ` flag specifies the BLAST algorithm to use. For the contamination filter, the goal is to identify and remove sequences with a very high similarity to the references. For this type of search it is best to use *megablast*, which is most appropriate for conducting within-species searches.

---------------
