'''
SuperCRUNCH: Rename_Merge module
                            
    Rename_Merge: The goal of this script is to replace the taxon names in a 
    sequence record description line with an updated/corrected name. If using the 
    Taxa_Assessment.py script, the fasta file to rename is called Unmatched_Taxa.fasta. 
    The user supplies a tab-delimited text file containing two columns, each with a binomial 
    name (genus and species separated by a space). Currently, this module does NOT support 
    the relabeling of subspecies names. The first column contains a taxon name that 
    should be corrected, and the second column contains the replacement name. The basis 
    of this replace file can be Unmatched_Records_Taxon_Names.log, which is an output from
    Taxa_Assessment.py. For each record, it is determined if the taxon name is in the 
    replacement file and, if so, the record will be updated with the new taxon name and 
    written to Renamed.fasta. Unless every name appearing in Unmatched_Records_Taxon_Names.log 
    is corrected, this new fasta file will inevitably be a subset of the records present 
    in Unmatched_Taxa.fasta.

    If the optional -m flag is used, the user must supply the full path to a fasta file to merge
    with the Renamed.fasta file. This should ideally be the Matched_Taxa.fasta produced by the
    Taxa_Assessment.py script. This merged fasta file should be used in the subsequent script,
    Parse_Loci.py.

[-r] The input taxon replacement file should be a tab-delimited text file with two columns
    of binomial names (genus and species separated by a space). The first column contains
    a taxon name that will be corrected, and the second column contains the replacement
    name. 
 
	Example of file structure for taxon replacement file:
	
    Varanus beccari	Varanus beccarii
    Varanus brevicauda	Varanus brevicaudus
    Varanus dumerili	Varanus dumerilii
    Varanus dumerilli	Varanus dumerilii
    Varanus eremias	Varanus eremius
    Vipera albizona	Montivipera albizona
    Vipera bornmuelleri	Montivipera bornmuelleri
    Vipera bulgardaghica	Montivipera bulgardaghica
    
[-m] An optional flag, should be the full path to the fasta file with valid taxon names.
    If using the Taxa_Assessment.py script, this would be called Matched_Taxa.fasta.
    This fasta file will be merged with the new fasta file containing renamed taxa, and
    is called 'Merged.fasta'. This should be used for the next script, Parse_Loci.py.

[-o] Full path to an existing directory in which the outputs will be written. 

	Outputs Files:
	
	Renamed.fasta - A fasta file containing all records with a matched taxon name. 

    Merged.fasta - A fasta file containing all records with a matched taxon name. 
    
    
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
import sys
import os
import argparse
import shutil
import time
import numpy as np
from Bio import SeqIO
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Rename_Merge: The goal of this script is to replace the taxon names in a 
    sequence record description line with an updated/corrected name. If using the 
    Taxa_Assessment.py script, the fasta file to rename is called Unmatched_Taxa.fasta. 
    The user supplies a tab-delimited text file containing two columns, each with a binomial 
    name (genus and species separated by a space). Currently, this module does NOT support 
    the relabeling of subspecies names. The first column contains a taxon name that 
    should be corrected, and the second column contains the replacement name. The basis 
    of this replace file can be Unmatched_Records_Taxon_Names.log, which is an output from
    Taxa_Assessment.py. For each record, it is determined if the taxon name is in the 
    replacement file and, if so, the record will be updated with the new taxon name and 
    written to Renamed.fasta. Unless every name appearing in Unmatched_Records_Taxon_Names.log 
    is corrected, this new fasta file will inevitably be a subset of the records present 
    in Unmatched_Taxa.fasta.

    If the optional -m flag is used, the user must supply the full path to a fasta file to merge
    with the Renamed.fasta file. This should ideally be the Matched_Taxa.fasta produced by the
    Taxa_Assessment.py script. This merged fasta file should be used in the subsequent script,
    Parse_Loci.py.
    
    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to a fasta file with taxon "
                            "names to replace; 'Unmatched_Taxa.fasta'")
    
    parser.add_argument("-r", "--replace",
                            required=True,
                            help="REQUIRED: The full path to a text file containing "
                            "all taxon names to be replaced, and the replacement names.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-m", "--merge",
                            required=False,
                            default=None,
                            help="OPTIONAL: The full path to a fasta file containing "
                            "valid taxon names; 'Matched_Taxa.fasta'")

    return parser.parse_args()

def index_fasta(f):
    """
    Use SeqIO index function to load fasta file (f), which
    is the best option for massive files.
    """
    tb = datetime.now()
    print("\nIndexing fasta file {}, this could take some time...".format(f))
    
    record_index = SeqIO.index(f, "fasta")
    
    tf = datetime.now()
    te = tf - tb
    print("\tTotal time to index fasta file: {0} (H:M:S)\n\n".format(te))

    return record_index

def parse_current_taxa(f):
    """
    Retrieve original taxon names from first column in user supplied 
    taxon renaming file (f). Returns as a list with all names in uppercase.
    """
    print("\nGathering list of taxon names to replace from {}.".format(f))
    taxa = []
    
    with open(f, 'r') as fh:
        lines = fh.readlines()
        taxa = [l.split('\t')[0].upper() for l in lines]
        
    print("\tFound {} taxa.\n".format(len(taxa)))
    
    return taxa

def make_taxa_dictionary(f):
    """
    Retrieve original and replacement taxon names from the two 
    columns in user supplied taxon renaming file (f). Returns 
    a dictionary where original name is key and replacement name
    is value.
    """
    taxa_dict = {}
    
    with open(f, 'r') as fh:
        for line in fh:
            parts = [l.strip().upper() for l in line.split('\t')]
            taxa_dict[parts[0]] = parts[1]
            
    return taxa_dict

def get_accession(line):
    """
    Retrieve the accession number from '>' line in a fasta file.
    """
    acc = [l.strip('>') for l in line.split()][0]
    
    return acc

def get_taxon(line):
    """
    Retrieve the taxon name from '>' line in a fasta file.
    Will fetch the species (binomial) name. Input line is 
    in uppercase so returned name is also uppercase.
    """
    if line.split() >= int(3):
        parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()][1:3]
        taxon = " ".join(parts)
        
    return taxon

def search_fasta_and_write(f, taxa_list, taxa_dict, record_index, out_f):
    """
    Input fasta file (f) is searched using a list of names recovered from the
    taxon renaming file (taxa_list), if name match is found then it is 
    replaced with an updated name from the taxon renaming file using a
    dictionary built (taxa_dict). To simplify, the original fasta file
    was indexed as a dictionary (record_index) and the output file (out_f) is 
    manually created from the new composite description line and the 
    original sequence (found by the accession number). Tracks progress.
    """
    tb = datetime.now()
    replace_cnt = int(0)
    record_cnt = int(0)
    
    with open(f, 'r') as fh_in:
        print("Replacing names and writing to new fasta file...")
        lcnt = int(0)
        for line in fh_in:
            lcnt += 1
            if lcnt % 10000 == 0:
                print("\tLine {}...".format(lcnt))

            if line.startswith('>'):
                if '|' in line:
                    line.replace("|"," ")
                record_cnt += 1 
                uline = line.upper().strip()
                taxon = get_taxon(uline)
                
                if taxon in taxa_list:
                    acc = get_accession(line)
                    parts = [l.strip() for l in line.split()]
                    assemble = [parts[0], taxa_dict[taxon].capitalize()]
                    for p in parts[3:]:
                        assemble.append(p)
                    newline = " ".join(assemble)
                    seq = record_index[acc].seq
                    
                    #Non-interleaved fasta is better for downstream scripts
                    #remember newline already contains > and full description
                    with open(out_f,'a') as fh_out:
                        fh_out.write("{}\n{}\n".format(newline, seq))
                    replace_cnt += 1
                    
    tf = datetime.now()
    te = tf - tb
    print("Relabeled {0} of {1} records.".format(replace_cnt, record_cnt))
    print("Total time to write fasta file with updated names: {0} (H:M:S)\n\n".format(te))

def find_upper(count):
    """
    Find where a number fits in a distribution. Will find 
    when the number is exceeded by a value in the distribution
    and returns that value. Needed for tracking progress on
    big lists (file lines, accession numbers, etc).
    """
    vala = np.array([100, 1000, 10000])
    valb = np.arange(100000, 1000000000, 100000)
    vals = np.append(vala,valb)
    for i in vals:
        if count <= i:
            upper = i
            break
        
    return upper

def merge(accs1, accs2, outname):
    """
    As it turns out writing a non-interleaved fasta file
    provides a major speed up in file searching for subsequent
    scripts, so merge fastas in this format.
    """
    acnt = int(0)
    upper1 = find_upper(len(accs1))
    upper2 = find_upper(len(accs2))
    tb = datetime.now()
    
    with open(outname, 'a') as fh_out:
        for a in accs1:
            fh_out.write(">{0}\n{1}\n".format(accs1[a].description, accs1[a].seq))
            acnt += 1
            interval = upper1 // 50
            if acnt % interval == 0:
                print("\t{} records written...".format(acnt))
                
        for a in accs2:
            fh_out.write(">{0}\n{1}\n".format(accs2[a].description, accs2[a].seq))
            acnt += 1
            interval = upper2 // 50
            if acnt % interval == 0:
                print("\t{} records written...".format(acnt))
                
    tf = datetime.now()
    te = tf - tb
    print("Wrote {0} records to merged fasta file".format(acnt))
    print("Total time to write merged fasta file: {0} (H:M:S)\n\n".format(te))
                    
def main():
    args = get_args()
    tb = datetime.now()
    os.chdir(args.outdir)
    
    taxa_list = parse_current_taxa(args.replace)
    taxa_dict = make_taxa_dictionary(args.replace)
    recs_rename = index_fasta(args.input)
    
    search_fasta_and_write(args.input, taxa_list, taxa_dict, recs_rename, "Renamed.fasta")
    
    if args.merge is not None:
        print("\nPreparing to merge fasta files.\n")
        recs_renamed = index_fasta("Renamed.fasta")
        recs_matched = index_fasta(args.merge)
        print("Writing merged fasta file...")
        merge(recs_renamed, recs_matched, "Merged.fasta")
        
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished renaming. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
                 
if __name__ == '__main__':
    main()
