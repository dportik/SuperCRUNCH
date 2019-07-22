'''
SuperCRUNCH: Infer_Supermatrix_Combinations module

    Infer_Supermatrix_Combinations: A tool for calculating how many supermatrix 
    combinations are available, given the number of filtered sequences available 
    for each taxon for each locus. If all taxa have only one sequence available, 
    the answer is one, but if taxa have multiple sequences available, this number
    will be extremely large. Relies on the [locus]_species_log.txt files produced 
    from the Filter_Seqs_and_Species.py module to calculate the number of 
    sequences available per taxon. The log files for all loci should be present
    in the input directory for the calculation to be accurate. No output files 
    are created, rather the information is logged to the screen.
    
    DEPENDENCIES: None.

-------------------------
Compatible with Python 2.7 & 3.7
Dependencies: 
    None
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
January 2019
Distributed under the 
GNU General Public Lincense
'''

import argparse
import os
from decimal import Decimal

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Infer_Supermatrix_Combinations: A tool for calculating how many supermatrix 
    combinations are available, given the number of filtered sequences available 
    for each taxon for each locus. If all taxa have only one sequence available, 
    the answer is one, but if taxa have multiple sequences available, this number
    will be extremely large. Relies on the [locus]_species_log.txt files produced 
    from the Filter_Seqs_and_Species.py module to calculate the number of 
    sequences available per taxon. The log files for all loci should be present
    in the input directory for the calculation to be accurate. No output files 
    are created, rather the information is logged to the screen.
    DEPENDENCIES: None.
	-----------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--in_dir",
                            required=True,
                            help="REQUIRED: The full path to a directory which "
                            "contains all the [locus]_species_log.txt files.")
    
    return parser.parse_args()

def parse_log(f):
    """
    Input file is tab delimited, six columns:
    Taxon Acession SeqLength Vouchered PassedTranslation SeqsAvailable
    Gather all entries for columns five and 1 and return as 
    list and set, respectively.
    """
    print("\tParsing information in {}".format(f))
    
    with open(f, 'r') as fh:
        #skip first line which is column labels
        next(fh)
        #if line not blank, split by tab, remove whitespace, take index 5 (col 6)
        seqs = [int((line.split('\t')[5].strip())) for line in fh if line.strip()]
        
    with open(f, 'r') as fh:
        #skip first line which is column labels
        next(fh)
        #if line not blank, split by tab, remove whitespace, take index 0 (col 1)
        taxa = set([(line.split('\t')[0].strip()) for line in fh if line.strip()])
        
    return seqs, taxa

def make_product(flist):
    """
    For every file in flist, obtains a list of number of sequences
    available and the set of taxon names. Adds these to the larger
    list of sequences and set of taxa. Multiplies all numbers in the
    sequence list to obtain the number of possible combinations.
    Returns combinations (integer), list of sequence numbers, and
    set of taxon names.
    """
    total_seqs = []
    taxa_set = set()
    
    for f in flist:
        seqs, taxa = parse_log(f)
        total_seqs.extend(seqs)
        taxa_set.update(taxa)
        
    seqs = sum(total_seqs)
    
    combinations = int(1)
    for n in total_seqs:
        combinations *= n
    return combinations, seqs, taxa_set

def main():
    args = get_args()
    os.chdir(args.in_dir)
    
    flist = sorted([f for f in os.listdir('.') if f.endswith('_species_log.txt')])
    print("\n\nFound {} loci to examine.\n\n".format(len(flist)))
    
    combinations, seqs, taxa = make_product(flist)
    
    print("\n\nFound {:,} total sequences for {:,} taxa.".format(seqs, len(taxa)))
    
    #Use special symbol to separate large number with commas.
    print("\n\nNumber of possible supermatrix combinations (unwieldy integer) = {:,}.\n"
              .format(combinations))
    
    #scientific notation, two formats
    val = '*10^'.join(((str("{:.2E}".format(Decimal(combinations)))).split('E+')))
    print("\nNumber of possible supermatrix combinations (scientific notation) = {:.2E}, or {}.\n\n"
              .format(Decimal(combinations), val))

if __name__ == '__main__':
    main()

