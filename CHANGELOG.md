## SuperCRUNCH Version History

- Version 1.2:
  - Made all modules compatible with Python 2.7 and Python 3.7.
  - SQL now implemented in `Parse_Loci.py` (30x speedup), `Filter_Seqs_and_Species.py` (3x speedup), and `Taxon_Assessment.py` (3x speedup).
  - Added output directory specification to all modules.
  - Two trimming modules now included: `Trim_Alignments_Trimal.py` and `Trim_Alignments_Custom.py`. The `Trim_Alignments_Custom.py` module allows finding start and stop block positions, and row-wise (internal) sliding window trimming based on divergence.
  - Added new module `Filter_Fasta_by_Min_Seqs.py` to filter fasta files using a minimum number of sequences.
  - Output directory structures improved for all modules.
  - Added `--quiet` option to `Filter_Seqs_and_Species.py` for less output on screen (useful when processing large numbers of loci).
  - Added option `--numerical` to `Fasta_Get_Taxa.py` to allow non-alphabetical identifiers for subspecies/trinomial name combinations. This allows museum, field, or numerical codes to be discovered.
  - Re-ordered tasks in `Cluster_Blast_Extract.py` to allow completion of all steps for one fasta file before moving to next fasta file in sequence.
  - Added multithreading for BLAST searches and new --bp_bridge flag for coordinate merging in `Cluster_Blast_Extract.py` and `Reference_Blast_Extract.py`.
  - Remove empty fasta files sometimes produced by `Coding_Translation_Tests.py`.
  - Complete code re-write for `Align.py`, `Cluster_Blast_Extract.py`, `Filter_Seqs_and_Species.py`, `Parse_Loci.py`, `Taxon_Assessment.py`.
  - Module `Relabel_Fasta.py` is now `Fasta_Relabel_Seqs.py`.
    
- Version 1.1:
  - Added multithreading option for MAFFT and Clustal-O in `Align.py`
  - Added multithreading option for MAFFT in `Adjust_Direction.py`
  - Added arg to specify output directory for `Concatenation.py`
  - Corrected output column labeling in label key output files from `Relabel_Fasta.py`
  - Added gappyout option for trimming with trimAl in `Trim_Alignments.py`
  - Output sequences failing similarity searches to own file in `Cluster_Blast_Extract.py` and `Reference_Blast_Extract.py`
  - Updated documentation on wiki pages
  