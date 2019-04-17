# Sequence Alignment

+ [Overview](#SA)
+ [Adjust_Direction.py](#AD)
+ [Coding_Translation_Tests.py](#CTT)
+ [Align.py](#A)

---------------

## **Overview** <a name="SA"></a>

![F4](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig4.jpg)

**SuperCRUNCH** includes two pre-alignment steps and several options for multiple sequence alignment. One pre-alignment step (`Adjust_Direction.py`) adjusts the direction of sequences and produces unaligned fasta files with all sequences written in the correct orientation. This step is necessary to avoid major pitfalls with aligners. Loci can be aligned using `Align.py` with **MAFFT**, **MUSCLE**, **Clustal-O**, or all these aligners sequentially. For coding loci, the **MACSE** translation aligner is also available, which is capable of aligning coding sequences with respect to their translation while allowing for multiple frameshifts or stop codons. To use this aligner the `Coding_Translation_Tests.py` module should be used to identify the correct reading frame of sequences, adjust them to the first codon position, and ensure completion of the final codon. Although MACSE can be run on a single set of reliable sequences (e.g., only those that passed translation), it has an additional feature allowing the simultaneous alignment of a set of reliable sequences and a set of unreliable sequences (e.g., those that failed translation) using different parameters. The `Coding_Translation_Tests.py` module can be used to generate all the necessary input files to perform this type of simultaneous alignment using **MACSE**.


---------------

### Adjust_Direction.py <a name="AD"></a>

The purpose of `Adjust_Direction.py` is to check sequences to ensure their proper direction before performing alignments with `Align.py` or additional filtering using `Coding_Translation_Tests.py`. `Adjust_Direction.py` is designed to work for a directory of fasta files, and it will adjust sequences in all unaligned fasta files and output unaligned fasta files with all sequences in the 'correct' direction. `Adjust_Direction.py` uses **MAFFT** to perform adjustments, and the default setting uses the *--adjustdirection* implementation of mafft. If the optional `--accurate` flag is included, it will use the *--adjustdirectionaccurately* option, which is slower but more effective with divergent sequences. 

The standard output file from **MAFFT** is an interleaved fasta file containing aligned sequences written in lowercase. In this file, sequences that have been reversed are flagged by writing an `_R_` at the beginning of the record ID. `Adjust_Direction.py` takes this output file and converts it to a cleaner format. The alignment is stripped, sequences are re-written sequentially and in uppercase, and the `_R_` is removed. Sequences that were reversed are recorded in a locus-specific log file, and the total counts of reversed sequences across all loci is written to a general log file, described below.

The resulting unaligned fasta files can be used for `Coding_Translation_Tests.py` or for `Align.py`.

#### Basic Usage:

```
python Adjust_Direction.py -i <input directory>
```

#### Argument Explanations:

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the unaligned fasta files. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `--accurate`

> **Optional**: Use the --adjustdirectionaccurately implementation in MAFFT, rather than --adjustdirection. It is slower but more accurate, especially for divergent sequences.

#### Example Uses:

```
python Adjust_Direction.py -i /bin/Adjust/
```
> Above command will adjust all unaligned fasta files in directory `Adjust/` using --adjustdirection in MAFFT.

```
python Adjust_Direction.py -i /bin/Adjust/ --accurate
```
> Above command will adjust all unaligned fasta files in directory `Adjust/` using --adjustdirectionaccurately in MAFFT.

Two outputs are created in the specified input directory for each fasta file, including:

+ `[fasta name]_Adjusted_Name_Log.txt`: Contains the full names of the sequences that were reversed in this particular fasta file, if any were adjusted.
+ `[fasta name]_Adjusted.fasta`: An unaligned fasta file which contains the correctly oriented sequences.

An additional output file is created:

+ `Log_Sequences_Adjusted.txt`: Contains the names of all fasta files and the number of sequences that were in the correct orientation or had to be reversed. Example contents:
```
Locus	Seqs_Correct_Direction	Seqs_Direction_Adjusted
12S	529	1
16S	565	13
CO1	484	0
CYTB	513	0
ND1	202	0
```

---------------

### Coding_Translation_Tests.py <a name="CTT"></a>

`Coding_Translation_Tests.py` can be used to identify translatable sequences from unaligned fasta files. This module performs the following tasks:
+ Sequences are translated in all forward frames to check for the presence of stop codons. 
    + If no stop codons are found for a forward frame, this frame is selected and the sequence passes translation.
    + If one stop codon is found and it is present in the final two codon positions of the sequence, this frame is selected and the sequence passes translation.
    + Sequences that have more than one stop codon in all frames fail the translation test. 
+ If a correct frame is identified, the sequence is adjusted so that the first base represents the first codon position. If no correct frame is found, the sequence is not adjusted.
+ For all sequences (pass and fail), the sequence length is adjusted with N's to ensure the final codon is complete. In other words, the total sequence length will be divisible by three.
+ Sequences are written as described above to corresponding output files. There are three output files written per input fasta file, described below.

If the `--rc` flag is included the translation will also be performed for the reverse complement, however if your sequences are all correctly oriented this is not recommended. The translation table should be specified with the `--table` flag. All NCBI translation table options are available, and can be selected using integers (*1-31*) or the shortcut terms provided (*standard, vertmtdna, invertmtdna, yeastmtdna, plastid*). If the `--table` flag is omitted, the default will be to use the Standard translation table. 

The `Coding_Translation_Tests.py` module can be used to determine if sequences are translatable, adjust sequences to the first codon position, and adjust sequence lengths to be divisible by three. Although these outputs were intended to be used for multiple sequence alignment with **MACSE**, the tests and outputs of `Coding_Translation_Tests.py` are likely to be useful for other purposes as well. 

#### Basic Usage:

```
python Coding_Translation_Tests.py -i <input directory> --table <translation table>
```

##### `-i <path-to-directory> `

> **Required**: The full path to a directory which contains the locus-specific fasta files to filter. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `--table <choice>`

> **Required**: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid*, or any integer *1-31*. Table will be used for all files found in the input directory.

##### `--rc`

> **Optional**: In addition to forward frames, use reverse complement for translation tests. Not recommended if direction of sequences has already been adjusted.

#### Example Uses:

```
python Adjust_Direction.py -i /bin/Translate/ --table vertmtdna
```
> Above command will perform translation tests for each unaligned fasta file in the directory `Translate/` using the vertebrate mitochondrial code.


An output directory called `Output_Translation_Fasta_Files/` is created in the input directory. For each fasta file included, the following output files are created:

+ `[Fasta name]_All.fasta`: Contains all length and/or position adjusted sequences, pass and fail.
+ `[Fasta name]_Passed.fasta`: Contains all length and position adjusted sequences that passed translation.
+ `[Fasta name]_Failed.fasta`: Contains all length adjusted sequences that failed translation.

An additional output file is created:

+ `Log_Sequences_Filtered.txt`:  A summary log file which indicates how many sequences passed translation and failed translation for each fasta file processed. Example contents:
```
Locus	Seqs_Passed	Seqs_Failed
CO1	479	5
CYTB	507	6
ND1	192	10
ND2	1005	0
ND4	548	3
```

---------------

### Align.py <a name="A"></a>

`Align.py` can be used to perform multiple sequence alignment for a directory of unaligned fasta files using ***MAFFT***, ***MUSCLE***, ***Clustal-O***, and/or ***MACSE***. 

The aligner is specified using the `-a` flag:
+ `-a mafft`: Align using ***MAFFT***.
+ `-a muscle`: Align using ***MUSCLE***.
+ `-a clustalo`: Align using ***Clustal-O***.
+ `-a all`: Align using ***MAFFT***, ***MUSCLE***, and ***Clustal-O*** (sequentially).
+ `-a macse`: Translation align using ***MACSE***. Details are provided below.

**NOTE:** The aligners can be run simultaneously on the same directory of unaligned fasta files without interfering with one another. This can speed up the alignment process if multiple aligners are being used, and is an alternative to the `-a all` option.

***MAFFT***, ***MUSCLE***, and ***Clustal-O*** can be used for all loci (coding and non-coding). The usage of each aligner invokes default settings or auto selection of the alignment strategy. For example: `mafft --auto` and `clustalo --auto`. These settings should be useful for a majority of alignments. The `--accurate` flag can be used to change the settings for ***MAFFT*** and ***Clustal-O*** to the following:
+  For ***Clustal-O***, this de-selects the --auto option and enables full distance matrix for guide-tree calculation, full distance matrix for guide-tree calculation during iteration, and --iter=5, in which the guide-tree and HMM each undergo 5 iterations, rather than only one: `clustalo --full --full-iter --iter=8`.
+ For ***MAFFT***, this option changes the default from auto select to use the FFT-NS-i strategy: `mafft --retree 2 --maxiterate 1000`.

The improved accuracy using these settings comes at the cost of longer run times for each aligner. However, this may be desirable for divergent sequences or difficult alignments.

***MACSE*** should only be used for coding loci, as it is a translation aligner. The use of `-a macse` requires using additional arguments and other optional arguments are available. The `--mpath` flag must be used to supply the full path to the ***MACSE*** JAR file (ideally V2). The `--table` flag can be used to specify the translation table using a shortcut term or an integer (*standard, vertmtdna, invertmtdna, yeastmtdna, plastid, 1-6, 9-16, 21-23*). Note that unlike previous steps, not all tables are available in ***MACSE***. Unless specified, the default table used is the Standard code. The optional `--mem` flag can be used to assign an amount of memory (in GB). The `--accurate` can also be used for ***MACSE*** v2, which invokes `-local_realign_init 0.9 -local_realign_dec 0.9`. The default for these arguments is normally 0.5., and changing this will slow down optimizations but increase alignment accuracy (sometimes considerably). These search settings will more closely resemble ***MACSE*** v1, which had default values of 1.0 for both.

An additional feature of ***MACSE*** is to include a set of reliable sequences (for example those that passed translation) and a set of less reliable sequences that are suspected to contain errors, and align both simultaneously with different penalty parameters. To use this feature in ***MACSE*** the `--pass_fail` flag can be used. However, to work the sequence sets must be contained in two fasta files that follow this naming format:

+ `[prefix]_Passed.fasta`: Fasta file of reliable sequences.
+ `[prefix]_Failed.fasta`: Fasta file of unreliable sequences.

These outputs are produced by default in the `Coding_Translation_Tests.py` module. The prefix portion of the name should ideally be the abbreviation of the gene/locus. If one file is missing, the `--pass_fail` will not work. 

Because ***MACSE*** can deal with frameshifts and sequence errors, it will insert an `!` character at corrected bp locations in the output alignment. A cleaned fasta file is created after the alignment is completed, in which all instances of `!` are replaced by `N`.

Output files vary between aligners but will be moved to output directories created in the main input directory, with details below. 


#### Basic Usage:

```
python Align.py -i <input directory> -a <aligner> 
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the unaligned fasta files. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `-a <choice>`

> **Required**: Specify whether alignment is by mafft, macse, muscle, or clustalo. If macse, MUST provide flags --mpath and --table. Choices = *mafft, macse, muscle, clustalo, all*.

##### `--mpath <path-to-executable>`

> **Required** for `-a macse`: Full path to a macse jar file (ideally MACSE v2).

##### `--table <choice>`

> **Required** for `-a macse`: Specifies translation table. Choices = *standard, vertmtdna, invertmtdna, yeastmtdna, plastid, 1-6, 9-16, 21-23*.

##### `--mem <integer>`

> **Optional** for `-a macse`: An integer for how much memory to assign to macse (in GB). Default = 1.

##### `--pass_fail`

> **Optional** for `-a macse`: Specifies macse to perform dual alignment. Files in -i directory must follow labeling format: NAME_Passed.fasta, NAME_Failed.fasta.

##### `--accurate`

> **Optional**: Specifies more thorough search settings (for mafft, clustalo, or macse).

#### Example Uses:

```
python Align.py -i /bin/to_align/ -a mafft
```
> Above command will align all fasta files in the `to_align/` directory using mafft with the --auto strategy.

```
python Align.py -i /bin/to_align/ -a mafft --accurate
```
> Above command will align all fasta files in the `to_align/` directory using mafft with the FFT-NS-i strategy.

```
python Align.py -i /bin/to_align/ -a all --accurate
```
> Above command will align all fasta files in the `to_align/` directory using clustalo, muscle and mafft. The `--accurate` flag changes settings using clustalo and mafft only.

```
python Align.py -i /bin/coding_loci/ -a macse --mpath /bin/programs/macse_v2.03.jar --table vertmtdna --mem 10 
```
> Above command will align all fasta files in the `coding_loci/` directory using the macse translation aligner (jar file `macse_v2.03.jar`) with the vertebrate mtdna table and 10GB of memory.

```
python Align.py -i /bin/coding_loci/ -a macse --mpath /bin/programs/macse_v2.03.jar --table vertmtdna --mem 10 --pass_fail
```
> Above command will align all paired fasta files (`[prefix]_Passed.fasta`, `[prefix]_Failed.fasta`) in the `coding_loci/` directory using the macse translation aligner (jar file `macse_v2.03.jar`) with the vertebrate mtdna table and 10GB of memory.

Depending on the alignment option selected, one or more of the directories will be created with the following contents:

+ **Output_CLUSTALO_Alignments/**
    + For each unaligned input fasta file, this directory contains a corresponding output alignment file labeled `[fasta name]_CLUSTALO_Aligned.fasta`. Results from `-a clustalo` and `-a all`.
+ **Output_MAFFT_Alignments/**
    + For each unaligned input fasta file, this directory contains a corresponding output alignment file labeled `[fasta name]_MAFFT_Aligned.fasta`. Results from `-a mafft` and `-a all`.
+ **Output_MUSCLE_Alignments/**
    + For each unaligned input fasta file, this directory contains a corresponding output alignment file labeled `[fasta name]_MUSCLE_Aligned.fasta`. Results from `-a muscle` and `-a all`.
+ **Output_MACSE_Alignments/**
    + For each unaligned input fasta file, this directory contains corresponding output files labeled `[fasta name]_AA.fasta`, `[fasta name]_NT.fasta`, and `[fasta name]_NT_Cleaned.fasta`. Results only from `-a macse`.

I ***strongly*** recommend using multiple aligners and comparing the results. This is arguably the most important step in creating phylogenetic data sets, and obtaining quality alignments is critical before performing any subsequent analyses.

---------------
