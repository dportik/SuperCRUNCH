# Post-Alignment Tasks

+ [Overview](#FFT)
+ [Relabel_Fasta.py](#RF)
+ [Trim_Alignments.py](#TAS)
+ [Fasta_Convert.py](#FC)
+ [Concatenation.py](#C)

---------------

## **Overview** <a name="FFT"></a>


![F5](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig5.jpg)

After multiple sequence alignment is complete, the alignment files are ready for relabeling (`Relabel_Fasta.py`), trimming (`Trim_Alignments.py`), format conversion (`Fasta_Convert.py`), and concatenation (`Concatenation.py`). There are different options for relabeling, depending on whether population-level data are gathered or species-level data are gathered. The relabeling step is required for trimming, and trimming is always recommended. The relabeled and trimmed fasta files can be converted into nexus and phylip format, and for species-level data sets (one representative sequence per taxon) the alignments can be easily concatenated into a phylogenetic supermatrix.

---------------

### Relabel_Fasta.py <a name="RF"></a>

The goal of `Relabel_Fasta.py` is to relabel all sequence records in all fasta files contained in a directory. The fasta files can be aligned or unaligned, but should have the original sequence description lines present.

There are three relabeling strategies available, specified using the `-r` flag: 
+ `-r species`: A taxon name is constructed from the description line. This generally corresponds to the genus and species if records are labeled properly. Subspecies labels can also be included using the optional `-s` flag.
+ `-r accession`: The accession number will be used to label the sequence record.
+ `-r species_acc`. The record is labeled by taxon and accession number. Subspecies labels can also be included using the optional `-s` flag.

In all strategies, any spaces are replaced by underscores. If the optional `-s` flag is included with a text file containing subspecies (trinomial) names, then the taxon component of each of the above options will include the trinomial if present.

The most appropriate relabeling strategy will depend on the final goals. Relabeling with `-r species` is essential for concatenating alignments to create a phylogenetic supermatrix. Relabeling with `-r species_acc` is the best option for population-level data sets, in which each species is represented by multiple sequences. This way, taxon names are present and the trailing accession numbers distinguish sequences belonging to the same taxon.

For each file relabeled, a corresponding labeling key is produced. This tab-delimited file contains the accession number, taxon label, and description line for each record. This serves as a reference key, such that the relevant information can be tracked for all relabeled samples.

Examples of how each relabeling option works is shown in greater detail below.

#### Basic Usage:

```
python Relabel_Fasta.py -i <input directory> -r <relabel option>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the unaligned or aligned fasta files. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

##### `-r <choice>`

> **Required**: The strategy for relabeling sequence records. Choices = *species, accession, species_acc*.

##### `-s <path-to-file>`

> **Optional**: The full path to a text file containing all subspecies names to cross-reference in the fasta file.


#### Example Uses:

```
python Relabel_Fasta.py -i bin/aligns_to_relabel/ -r species -s bin/subspecies_list.txt
```
> Above command will relabel description lines using the taxon name, and will use the appropriate subspecies labels from `subspecies_list.txt` if they are present in the records. Action is performed for all fasta files located in the `aligns_to_relabel/` directory.

```
python Relabel_Fasta.py -i bin/aligns_to_relabel/ -r acc 
```
> Above command will relabel description lines using the accession number. Action is performed for all fasta files located in the `aligns_to_relabel/` directory.

```
python Relabel_Fasta.py -i bin/aligns_to_relabel/ -r species_acc -s bin/subspecies_list.txt
```
> Above command will relabel description lines using the taxon name plus accession number, and will use the appropriate subspecies labels from `subspecies_list.txt` if they are present in the records. Action is performed for all fasta files located in the `aligns_to_relabel/` directory.

#### Relabeling details:

Here I show how the relabeling works, in greater detail.

Let's use the following set of description lines from a fasta file:
```
>JN881132.1 Daboia russelii activity-dependent neuroprotector (ADNP) gene, partial cds
>KU765220.1 Sceloporus undulatus voucher ADL182 activity-dependent neuroprotector (adnp) gene, partial cds
>DQ001790.1 Callisaurus draconoides carmenensis isolate RWM 1480 cytochrome b (cytb) gene
```
Notice that the third entry has a subspecies label, but the first two entries do not.

Using `-r species` would produce:
```
>Daboia_russelii
>Sceloporus_undulatus
>Callisaurus_draconoides
```
The subspecies label in the third entry is ignored.

Using `-r accession` would produce:
```
>JN881132.1
>KU765220.1
>DQ001790.1
```
Using `-r species_acc` would produce:
```
>Daboia_russelii_JN881132.1
>Sceloporus_undulatus_KU765220.1
>Callisaurus_draconoides_DQ001790.1
```
Again, the subspecies label in the third entry is ignored.

If I supplied a subspecies text file using the `-s` flag that contained `Callisaurus draconoides carmenensis`, then the following would be produced with `-r species`:
```
>Daboia_russelii
>Sceloporus_undulatus
>Callisaurus_draconoides_carmenensis
```
The subspecies label in the third entry is now included.

And this would be produced using `-r species_acc`:
```
>Daboia_russelii_JN881132.1
>Sceloporus_undulatus_KU765220.1
>Callisaurus_draconoides_carmenensis_DQ001790.1
```
The subspecies label in the third entry is now included.

Depending on the relabeling strategy selected, one of the directories will be created with the following contents:

+ **Relabeled_Fasta_Files_Species/**
    + For each locus, this directory contains a corresponding fasta file labeled `[fasta name]_relabeled.fasta` and `[fasta name]_label_key.txt`. Results from `-r species`.
+ **Relabeled_Fasta_Files_Accession/**
    + For each locus, this directory contains a corresponding fasta file labeled `[fasta name]_relabeled.fasta` and `[fasta name]_label_key.txt`. Results from `-r accession`.
+ **Relabeled_Fasta_Files_SpeciesAccession/**
    + For each locus, this directory contains a corresponding fasta file labeled `[fasta name]_relabeled.fasta` and `[fasta name]_label_key.txt`. Results from `-r species_acc`.

Example contents from a `[fasta name]_label_key.txt` file:
```
Taxon	Accession	Description
KU097508.1	Acanthocercus adramitanus	Acanthocercus adramitanus isolate IBES10359 16S ribosomal RNA gene, partial sequence; mitochondrial
MG700133.1	Acanthocercus annectens	Acanthocercus annectens voucher USNM:589391 16S ribosomal RNA gene, partial sequence; mitochondrial
JX668132.1	Acanthocercus atricollis	Acanthocercus atricollis voucher EBG 2167 16S ribosomal RNA gene, partial sequence; mitochondrial
JX668138.1	Acanthocercus cyanogaster	Acanthocercus cyanogaster voucher MVZ 257937 16S ribosomal RNA gene, partial sequence; mitochondrial
JX668140.1	Acanthocercus yemensis	Acanthocercus yemensis voucher MVZ 236454 16S ribosomal RNA gene, partial sequence; mitochondrial
NC_014175.1	Acanthosaura armata	Acanthosaura armata mitochondrion, complete genome
MG935713.1	Acanthosaura crucigera	Acanthosaura crucigera voucher USNM:587019 16S ribosomal RNA gene, partial sequence; mitochondrial
```

---------------

### Trim_Alignments.py <a name="TAS"></a>

The alignments may require some amount of trimming to remove overhanging ends, poorly aligned regions, etc. `Trim_Alignments.py` simplifies this process by automating trimming for a directory of input files. `Trim_Alignments.py` relies on **trimAl**, and four options are available for trimming. Specifying `-a gt` will use the gap threshold method, and although the default value is 0.05 this can be changed using the `--gt ` flag to set a value between 0 and 1. The gt value is the minimum fraction of sequences without a gap required to retain an alignment column. Specifying `-a noallgaps` uses the noallgaps method, which simply removes any alignment column composed entirely of gaps. Specifying `-a both` runs the gap threshold method followed by the noallgaps method. Finally, specifying `-a gappyout` runs the gappyout method.

The input files can be in fasta, nexus, or phylip format, and the format is automatically detected by **trimAl**. The output format must be specified by the `-f ` flag, and includes fasta, nexus, or phylip format. Although any output format can be chosen, fasta format is recommended if the `Fasta_Convert.py` or `Concatenation.py` modules will be used.

**Note:** If the input files have not been relabeled at this point (ie they are aligned fasta files containing the original description lines), the identifier names appearing in the output files will default to accession numbers only. 

#### Basic Usage:

```
python Trim_Alignments.py -i <input directory> -f <output format> -a <trimal method>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the input alignment files. File formats can include fasta, nexus, or phylip, with one of the corresponding file extensions: '.fasta', '.fa', '.nexus', '.nex', '.phylip', or '.phy'. 

##### `-f <choice>`

> **Required**: Specifies the output file format for trimmed alignments. Choices = *fasta, nexus, phylip*.

##### `-a <choice>`

> **Required**: Specifies the trimal method for trimming alignments. Choices = *gt, noallgaps, both, gappyout*.

##### `--gt <value>`

> **Optional**: Specifies the gap threshold (gt) value for trimal, the minimum fraction of sequences without a gap. Must be between 0 and 1. Default = 0.05.


#### Example Uses:

```
python Trim_Alignments.py -i /bin/Trim/ -f fasta -a gt --gt 0.1
```
> Above command will trim all alignments present in the directory `Trim/` using the gap threshold method with a value of 0.1, and output files in fasta format.

```
python Trim_Alignments.py -i /bin/Trim/ -f phylip -a noallgaps
```
> Above command will trim alignments present in the directory `Trim/` using the noallgaps method and output files in phylip format.

```
python Trim_Alignments.py -i /bin/Trim/ -f nexus -a both --gt 0.08
```
> Above command will trim alignments present in the directory `Trim/` using the gap threshold method with a value of 0.08, and output files in nexus format.

```
python Trim_Alignments.py -i /bin/Trim/ -f fasta -a gappyout
```
> Above command will trim alignments present in the directory `Trim/` using the gappyout method and output files in fasta format.

---------------

### Fasta_Convert.py <a name="FC"></a>

The goal of `Fasta_Convert.py` is to convert a directory of aligned fasta files into both phylip and nexus formats. Note this should only be used after records have been renamed with species labels (format `genus_species` or `genus_species_subspecies`) or accession numbers, as the original NCBI description lines will cause severe issues with writing other file formats.


#### Basic Usage:

```
python Fasta_Convert.py -i <input directory>
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory which contains the ALIGNED fasta files. Fasta files in the directory must have extensions '.fasta' or '.fa' to be read.

#### Example Uses:

```
python Fasta_Convert.py -i /bin/relabeled_alignments/
```
> Above command will convert all fasta alignments in the directory `relabeled_alignments/` into nexus and phylip format.

Two output folders are created in the input directory containing the fasta files. The directory labels and their contents are described below:

+ **Output_Phylip_Files/**
    + For each locus, this directory contains a corresponding phylip file, labeled `[fasta name].phy`.
+ **Output_Nexus_Files/**
    + For each locus, this directory contains a corresponding nexus file, labeled `[fasta name].nex`.


---------------

### Concatenation.py <a name="C"></a>

The goal of using `Concatenation.py` is to combine many alignments into a single concatenated alignment. 

To work properly, sequence names should have been relabeled across all files and should ideally take the format of `genus_species` or `genus_species_subspecies`, in other words taxon names. However, any sequence label can be used as long as the labels are consistent across alignments. There cannot be any duplicate sequence labels within any of the alignment files, or an error will be thrown: `ValueError: Duplicate key [sequence label]`. The complete set of taxa/labels is inferred from all the alignments, for each taxon/label the sequences from all loci are retrieved. If a taxon/label is absent from an alignment, a missing sequence is generated using the symbol selected with the `-s` flag (options: N, dash, or ? symbol). The input files can be in fasta or phylip format, specified using the `--in_format` flag. The output format must be specified as fasta or phylip using the `--out_format` flag. Taxa/labels are written in alphabetical order in the concatenated alignment. The arrangement of sequences in the concatenated alignment corresponds to the alphabetical order of the names of the input alignment files. Besides producing the concatenated alignment, two additional outputs are produced and described below. The outputs are written to a directory specified by the user with the `-o` flag.
 
This script takes advantage of python dictionary structures to read alignments and store sequences. As a result, it works extremely fast and is more than capable of concatenating thousands of large alignment files. 


#### Basic Usage:

```
python Concatenation.py -i <input directory> -o <output dirctory> --in_format <input alignment format> --out_format <output alignment format> -s <missing data symbol> 
```

##### `-i <path-to-directory>`

> **Required**: The full path to a directory containing the aligned files. 

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write the concatenated alignment and other output files.

##### `--in_format <choice>`

> **Required**: The file format of the INPUT alignments. Choices = *fasta, phylip*. The files must have one of the corresponding extensions to be read: '.fasta', '.fa', '.phylip', or '.phy'.

##### `--out_format <choice>`

> **Required**: The file format for the OUTPUT concatenated alignment. Choices = *fasta, phylip, nexus*.

##### `-s <choice>`

> **Required**: A base pair symbol used to represent missing data sequences. Choices = *dash, N, or ?*

#### Example Uses:

```
python Concatenation.py -i /bin/Final_Alignments/ -o /bin/Concatenated/ --in_format fasta --out_format phylip -s dash 
```
> Above command will concatenate all the fasta alignments in the directory `Final_Alignments/` and produce a phylip file in `Concatenated/`. Missing data are represented with the `-` symbol.

```
python Concatenation.py -i /bin/Final_Alignments/ -o /bin/Concatenated/ --in_format phylip --out_format fasta -s ?
```
> Above command will concatenate all the phylip alignments in the directory `Final_Alignments/` and produce a fasta file in `Concatenated/`. Missing data are represented with the `?` symbol.

Three output files are created in the specified output directory, including:

+ `Concatenated_Alignment.fasta`,`Concatenated_Alignment.phy`, or `Concatenated_Alignment.nex`: The final concatenated alignment in the output format specified.
+ `Data_Partitions.txt`: Text file displaying the order in which loci were concatenated and their corresponding base pairs within the alignment. These are the data partitions or charsets, and can be used for model testing or partitioning in analyses.
+ `Taxa_Loci_Count.log`: A count of the number of sequences that were available for each taxon/label, in other words the number of alignments the taxon/label was present in. If taxon names were used, this will be a full list of the species included in the concatenated alignment.


---------------