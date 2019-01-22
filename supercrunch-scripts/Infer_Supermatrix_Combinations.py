import argparse
import os
from decimal import Decimal

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Infer_Supermatrix_Combinations: A tool for calculating how many supermatrix combinations are available,
    given the number of filtered sequences available for each taxon for each locus. If all taxa have only
    one sequence available, the answer is one, but if taxa have multiple sequences available, this number
    will be extremely large. Relies on the [locus]_species_log.txt files produced from the Filter_Seqs_and_Species.py
    module to calculate the number of sequences available per taxon. The log files for all loci should be present
    in the input directory for the calculation to be accurate. No output files are created, rather the information is logged to the screen.
    		
    DEPENDENCIES: None.
	-----------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains all the [locus]_species_log.txt files.")
    return parser.parse_args()

def parse_log(f):
    '''
    Input file is tab delimited, five columns - Taxon	Accession	SeqLength	PassedTranslation	SeqsAvailable
    Gather all entries for column five (SeqsAvailable) and return list.
    '''
    print "\tParsing information in {}".format(f)
    with open(f, 'r') as fh:
        #skip first line which is column labels
        next(fh)
        #if line not blank, split by tab, remove whitespace, take index 4
        seqs = [int((line.split('\t')[4].strip())) for line in fh if line.strip()]
    with open(f, 'r') as fh:
        #skip first line which is column labels
        next(fh)
        #if line not blank, split by tab, remove whitespace, take index 4
        taxa = [(line.split('\t')[0].strip()) for line in fh if line.strip()]
    return seqs, taxa

def make_product(f_list):
    total_seqs = []
    total_taxa = []
    
    for f in f_list:
        seqs, taxa = parse_log(f)
        total_seqs.extend(seqs)
        total_taxa.extend(taxa)
        
    taxa = set(total_taxa)
    
    seqs = sum(total_seqs)
    
    combinations = int(1)
    for n in total_seqs:
        combinations *= n
    return combinations, seqs, taxa

def main():
    args = get_args()
    os.chdir(args.in_dir)
    f_list = sorted([f for f in os.listdir('.') if f.endswith('_species_log.txt')])
    print "\n\nFound {} loci to examine.\n\n".format(len(f_list))
    combinations, seqs, taxa = make_product(f_list)
    print "\n\nFound {:,} total sequences for {:,} taxa.".format(seqs, len(taxa))
    #Use special symbol to separate large number with commas.
    print "\n\nNumber of possible supermatrix combinations (unwieldy integer) = {:,}.\n".format(combinations)
    #scientific notation, two formats
    val = '*10^'.join(((str("{:.2E}".format(Decimal(combinations)))).split('E+')))
    print "\nNumber of possible supermatrix combinations (scientific notation) = {:.2E}, or {}.\n\n".format(Decimal(combinations), val)

if __name__ == '__main__':
    main()

