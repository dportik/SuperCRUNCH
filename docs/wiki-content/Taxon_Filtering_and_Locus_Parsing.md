# Taxon Filtering and Locus Parsing

+ [Overview](#TFLE)
+ [Taxa_Assessment.py](#TA)
+ [Rename_Merge.py](#RM)
+ [Parse_Loci.py](#PL)

---------------

## **Overview** <a name="TFLE"></a>

![F1](https://github.com/dportik/SuperCRUNCH/blob/master/docs/Fig1.jpg)

This section deals with the first filtering steps, which include screening sequence records for taxon names and locus identity. To run these steps in **SuperCRUNCH**, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms. Information on obtaining these files is provided in the section above. The `Taxa_Assessment.py` and `Rename_Merge.py` scripts are optional, but highly recommended. `Taxa_Assessment.py` identifies all valid and invalid taxon names contained within the starting fasta file. `Rename_Merge.py` is an optional data cleaning step that can be used to replace invalid taxon names with updated valid names for corresponding  sequence records. This relabeling step allows these records to pass the taxonomy filter in `Parse_Loci.py`, rather than be discarded. The original or updated fasta file is processed using `Parse_Loci.py`. For each locus included, a fasta file is produced. These fasta files contain records with valid taxon names that matched one or more of the search terms for the locus.

---------------

### Taxa_Assessment.py <a name="TA"></a>

The goal of this script is to search through records in a fasta file of NCBI nucleotide sequences to determine whether or not they contain a taxon name present in the description line that matches a taxon name in the user-supplied taxon list. The taxon names list can contain a mix of species (binomial name) and subspecies (trinomial name) labels.

Two output fasta files are written to the specified output directory: one containing only records with valid taxon names (`Matched_Taxa.fasta`), and one containing records with invalid taxon names (`Unmatched_Taxa.fasta`). The accession numbers for each of these fasta files are also written to separate files (`Matched_Records_Accession_Numbers.log`, `Unmatched_Records_Accession_Numbers.log`). Two log files are written which contain lists of the valid (`Matched_Taxon_Names.log`) and invalid taxon names (`Unmatched_Taxon_Names.log`) found across all records. The `Unmatched_Taxon_Names.log` file can be used to create the base file needed to relabel taxa in the `Rename_Merge.py` script. 

The decision to include or exclude subspecies labels is up to the user, and can be specified using the `--no_subspecies` flag. For a thorough explanation of how taxonomy searches are conducted and how this flag affects this step (and others), please see below. For all searches the sequence description line and the supplied taxon name are converted to uppercase, so the list of taxon names is not case-sensitive.

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
> Above command will search sequence records in `Start_Seqs.fasta` to find taxon names present in `Taxa_List.txt`, the output is written to the specified directory.

```
python Taxa_Assessment.py -i bin/Analysis/Start_Seqs.fasta -t bin/Analysis/Taxa_List.txt -o bin/Analysis/Output/ --no_subspecies
```
> Above command will search sequence records in `Start_Seqs.fasta` to find taxon names present in `Taxa_List.txt`, the output is written to the specified directory. This search will exclude the subspecies component of taxon labels in the taxon names file and fasta file.


#### Taxonomy searches with and without subspecies:

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

In the example above, `Leycesteria crocothyrsos` and `Linnaea borealis` would be the names used. Essentially, the `--no_subspecies` flag 'sinks' all the subspecies, and they are all lumped under the relevant species label.

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

The goal of this script is to search through records in a fasta file and replace invalid taxon names with new names that are compatible with the taxon list. The invalid taxon names file (`Unmatched_Taxon_Names.log`) from the `Taxa_Assessment.py` can be used to help create the replacement names file, as it contains all the taxon names that failed the taxon filter. A second column can be added, and as each name is inspected a replacement name can be added to the second column, or if the name cannot be rescued then the row can be deleted. The final replacement names file should be two-columns in tab-delimited format, and an example is shown below. Similar to other input text files, make sure to open the replacement names file in a text editor and ensure the format includes Unix line breaks (line breaks marked by `\n`, rather than `\r\n`) and UTF-8 encoding, otherwise extra characters may interfere with parsing the file correctly with **SuperCRUNCH**.

If the optional `-m` flag is used, records that are successfully relabeled are merged with another fasta file. Ideally, this fasta file should be the one containing all the records with valid taxon names (*Matched_Taxa.fasta*). The resulting merged fasta file (*Merged.fasta*) should then be used for the `Parse_Loci.py` script.

#### Basic Usage:

```
python Rename_Merge.py -i <fasta file> -r <taxon renaming file> -o <output directory>
```

#### Argument Explanations:

##### `-i <path-to-file>`

> **Required**: The full path to a fasta file with taxon names to replace (`Unmatched_Taxa.fasta`).

##### `-r <path-to-file>`

> **Required**: The full path to a two-column text file containing all taxon names to be replaced, and the replacement names.

##### `-o <path-to-directory>`

> **Required**: The full path to an existing directory to write output files.

##### `-m <full-path-to-file>`

> **Optional**: The full path to a fasta file containing valid taxon names (`Matched_Taxa.fasta`). 

#### Example Use:

```
python Rename_Merge.py -i bin/Rename/Unmatched_Taxa.fasta -r bin/Rename/taxon_relabeling.txt -o bin/Rename/Output/
```
> Above command will attempt to rename sequence records in `Unmatched_Taxa.fasta` following the file `taxon_relabeling.txt`, and the output is written to the specified directory.
```
python Rename_Merge.py -i bin/Rename/Unmatched_Taxa.fasta -r bin/Rename/taxon_relabeling.txt -o bin/Rename/Output/ -m bin/Rename/Matched_Taxa.fasta
```
> Above command will attempt to rename sequence records in `Unmatched_Taxa.fasta` following the file `taxon_relabeling.txt`, and merge these relabeled records with those in `Matched_Taxa.fasta`. The output is written to the specified directory.


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

In other cases, taxon names may have been updated and now represent synonymies, or may have been accidentally misspelled. Using a organism-specific taxonomy browser can help clarify these situations. Synonymies, misspellings, and name changes represent examples of records that are worth rescuing through relabeling, and using `Rename_Merge.py` to do so will result in higher quality data. 

---------------

### Parse_Loci.py <a name="PL"></a>

The goal of `Parse_Loci.py` is to search through the fasta file of starting sequences and identify records for each gene/marker/locus. Each record that matches a locus must also contain a valid taxon name to be retained. To run `Parse_Loci.py`, you will need to provide a fasta file of downloaded nucleotide sequence records, a file containing a list of taxa, and a file containing a list of loci and associated search terms.  

The taxon names list can contain a mix of species (binomial) and subspecies (trinomial) names. Detailed instructions for the format of this file is provided in the `Taxa_Assessment.py` section. The optional `--no_subspecies` flag can be used, and its effect is also described in great detail in the `Taxa_Assessment.py` section.

The locus file contains the search terms that are used to identify matching records. Detailed information about the format of this file can be found in the **Obtaining Loci and Search Terms** section. 
 
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
python Parse_Loci.py -i bin/Loci/Merged.fasta -l bin/Loci/locus_search_terms.txt -t bin/Loci/Taxa_List.txt -o bin/Loci/Output/
```
> Above command will use the locus search terms file `locus_search_terms.txt` and the taxon names file `Taxa_List.txt` to parse records in `Merged.fasta`, writing outputs to the specified directory.
```
python Parse_Loci.py -i bin/Loci/Merged.fasta -l bin/Loci/locus_search_terms.txt -t bin/Loci/Taxa_List.txt -o bin/Loci/Output/ --no_subspecies
```
> Above command will use the locus search terms file `locus_search_terms.txt` and the taxon names file `Taxa_List.txt` to parse records in `Merged.fasta`, ignoring the subspecies component of taxon names. Outputs are written to the specified directory.


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

Identifying loci in the sequence records through word matching is not expected to be perfect, and there may be non-target sequences in the resulting fasta files. For this reason, the fasta files should be subjected to orthology filtering before attempting to create sequence alignments or performing analyses, as described in the next section.

---------------

