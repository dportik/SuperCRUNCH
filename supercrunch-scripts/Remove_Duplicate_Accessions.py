'''
Usage: python Remove_Duplicate_Accessions.py -i [fasta file] REQUIRED
							                 -o [path to output directory] REQUIRED
                            
    Remove_Duplicate_Accessions: Search through a GenBank fasta file and remove any duplicate records (based
    on identical accession numbers). This prevents errors from happening in BioPython when a fasta
    file is attempted to be read in as a dictionary, as no duplicate keys are allowed. The generator
    function used here makes it efficient to run this for massively large files. 

-------------------------
Written for Python 2.7
Python modules required:
	-BioPython (using SeqIO module)
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated November 2018
'''
import argparse
import os
from Bio import SeqIO

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""------------------------------------------------------------------------------
    Remove_Duplicate_Accessions: Search through a GenBank fasta file and remove any duplicate records (based
    on identical accession numbers). This prevents errors from happening in BioPython when a fasta
    file is attempted to be read in as a dictionary, as no duplicate keys are allowed. The generator
    function used here makes it efficient to run this for massively large files. 
    DEPENDENCIES: Python: BioPython.
	------------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input", required=True, help="REQUIRED: The full path to a fasta file of GenBank sequence data to filter")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: The full path to an existing directory to write output fasta file.")
    return parser.parse_args()

def remove_dup_seqs(records):
    print "\n\nSearching through fasta file..."
    accs = set()
    for record in records:
        a = record.id
        if a in accs:
            print "Ignoring {}".format(record.id)
            continue
        accs.add(a)
        yield record

def main():
    args = get_args()
    os.chdir(args.out_dir)
    records = remove_dup_seqs(SeqIO.parse(args.input, "fasta"))
    write_fasta = SeqIO.write(records, "Cleaned.fasta", "fasta")
    
if __name__ == '__main__':
    main()
