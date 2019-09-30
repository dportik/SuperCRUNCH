'''
SuperCRUNCH: Remove_Duplicate_Accessions module
                            
    Remove_Duplicate_Accessions: Search through a fasta file and remove any duplicate 
    records, based on identical accession numbers. This prevents errors from happening 
    in BioPython when a fasta file is attempted to be read in as a dictionary, as no 
    duplicate keys are allowed. The generator function used here makes it efficient 
    to run this for massively large files. 

-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
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
from Bio import SeqIO

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""------------------------------------------------------------------------------
    Remove_Duplicate_Accessions: Search through a fasta file and remove any duplicate 
    records, based on identical accession numbers. This prevents errors from happening 
    in BioPython when a fasta file is attempted to be read in as a dictionary, as no 
    duplicate keys are allowed. The generator function used here makes it efficient 
    to run this for massively large files. 

    DEPENDENCIES: Python: BioPython.
	------------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to a fasta file of "
                            "sequence data to filter.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output fasta file.")
    
    return parser.parse_args()

def remove_dup_seqs(records):
    """
    A generator that builds a set of accession numbers
    and returns records that are not already found in
    the accession set. 
    """
    print("\n\nSearching through fasta file...")
    accs = set()
    for record in records:
        a = record.id
        if a in accs:
            print("Ignoring {}".format(record.id))
            continue
        accs.add(a)
        yield record

def main():
    args = get_args()
    os.chdir(args.outdir)
    prefix = args.input.split('/')[-1].split('.')[0]
    records = remove_dup_seqs(SeqIO.parse(args.input, "fasta"))
    write_fasta = SeqIO.write(records, "{}-Cleaned.fasta".format(prefix), "fasta")
    
if __name__ == '__main__':
    main()
