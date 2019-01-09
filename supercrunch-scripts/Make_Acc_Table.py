'''
Usage: python Make_Acc_Table.py -i [directory containing fasta files] REQUIRED
                                -s [full path to text file with subspecies names to include] OPTIONAL


    Make_Acc_Table: A tool for extracting taxon and accession number information from
    a directory of fasta files containing records with the original description lines.
    For example:
    >FJ745388.1 Abelia chinensis voucher N. Pyck 1989-2220 internal transcribed spacer 1, partial sequence;

    A comprehensive list of taxa is collected across all fasta files in the directory. The taxon
    names to search are binomial by default (only genus + species), but if the optional -s flag 
    is used then the subspecies names contained in the supplied text file will also be searched for 
    and included (if present). The table is written with species in the first column (sorted 
    alphabetically) and the fasta files composing the remaining columns (labeled by fasta name), 
    with rows filled with accession numbers or dashes depending on the presence or absence of the 
    taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files. 

    To be recognized, a fasta file must be labeled as "NAME.fasta" or "NAME.fa". The NAME portion 
    should not contain any periods or spaces, but can contain underscores. 

-------------------------
Written for Python 2.7
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated December 2018
'''
import os
import shutil
import argparse
import random
from datetime import datetime

def get_args():
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Make_Acc_Table: A tool for extracting taxon and accession number information from
    a directory of fasta files containing records with the original description lines.
    For example:
    '>FJ745388.1 Abelia chinensis voucher N. Pyck 1989-2220 internal transcribed spacer 1, partial sequence;'.

    A comprehensive list of taxa is collected across all fasta files in the directory. The taxon
    names to search are binomial by default (only genus + species), but if the optional -s flag 
    is used then the subspecies names contained in the supplied text file will also be searched for 
    and included (if present). The table is written with species in the first column (sorted 
    alphabetically) and the fasta files composing the remaining columns (labeled by fasta name), 
    with rows filled with accession numbers or dashes depending on the presence or absence of the 
    taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files.

    To be recognized, a fasta file must be labeled as "NAME.fasta" or "NAME.fa". The NAME portion 
    should not contain any periods or spaces, but can contain underscores. 
    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files.")
    parser.add_argument("-s", "--subspecies", required=False, default=None, help="OPTIONAL: The full path to a text file containing all subspecies names to cross-reference in the fasta file.")
    return parser.parse_args()        

def parse_taxa(f):
    print "\nParsing taxon information from {}.".format(f)
    with open(f, 'r') as fh_f:
        subspecies_set = set([line.upper().strip().replace(" ","_") for line in fh_f if len(line.split()) == int(3)])
    subspecies = sorted(subspecies_set)
    print "\tFound {} subspecies names to include.\n".format(len(subspecies))
    return subspecies

def get_taxon(line):
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split() if line.split() >= int(3)][1:3]
    taxon_sp = "_".join(parts1)
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split() if line.split() >= int(4)][1:4]
    taxon_ssp = "_".join(parts2)
    return taxon_sp, taxon_ssp

def get_taxa_fasta(f, subspecies):
    '''
    Retrieve the names following the > in all lines
    from a fasta file, return as a list.
    '''
    taxa_set = set()
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith(">"):
                if '|' in line:
                    line.replace("|"," ")
                taxon_sp, taxon_ssp = get_taxon(line)
                if subspecies is not None:
                    if taxon_ssp in subspecies:
                        taxa_set.add(taxon_ssp)
                    else:
                        taxa_set.add(taxon_sp)
                else:
                    taxa_set.add(taxon_sp)
    taxa = sorted(taxa_set)
    return taxa
    
def collect_taxa(f_list, subspecies):
    '''
    Collect taxon names from all alignment files
    and return a sorted list of all unique names.
    '''
    taxon_set = set()
    for f in f_list:
        f_taxa = get_taxa_fasta(f, subspecies)
        taxon_set.update(set(f_taxa))
    taxon_list = sorted(taxon_set)
    return taxon_list

def f_taxon_dict(f, subspecies):
    '''
    Function to convert fasta file into
    dictionary structure with taxon as key
    and accession number as value.
    '''
    f_dict = {}
    with open(f, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    for line in lines:
        if line.startswith(">"):
            if '|' in line:
                line.replace("|"," ")
            acc = [l.strip('>') for l in line.split()][0]
            taxon_sp, taxon_ssp = get_taxon(line)
            if subspecies is not None:
                if taxon_ssp in subspecies:
                    f_dict[taxon_ssp] = acc
                else:
                    f_dict[taxon_sp] = acc
            else:
                f_dict[taxon_sp] = acc
    return f_dict

def taxa_acc_dict(dict_list, taxa):
    writing_dict = {}
    for taxon in taxa:
        writing_dict[taxon] = ""
        for d in dict_list:
            if taxon in d:
                writing_dict[taxon] += "{}\t".format(d[taxon])
            elif taxon not in d:
                writing_dict[taxon] += "{}\t".format("-")
    return writing_dict

def write_acc_table(taxa, writing_dict, f_names):
    with open("GenBank_Accession_Table.txt", 'a') as fh_out:
        fh_out.write('Taxon\t')
        for f in f_names:
            fh_out.write('{}\t'.format(f))
        fh_out.write('\n')
        for taxon in taxa:
            fh_out.write("{0}\t{1}\n".format(taxon.capitalize().replace("_"," "), writing_dict[taxon].strip('\t')))

def main():
    tb = datetime.now()
    args = get_args()
    subspecies = None
    if args.subspecies is not None:
        subspecies = parse_taxa(args.subspecies)
        
    os.chdir(args.in_dir)
    f_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])
    print "\nFound {0} files to use for accession table.".format(len(f_list))
    f_names = [f.split(".")[0] for f in f_list]
    
    taxa = collect_taxa(f_list, subspecies)
    print "\nFound {} unique taxa across fasta files.".format(len(taxa))
    print "\nGathering accession numbers for all taxa."
    dict_list = [f_taxon_dict(f, subspecies) for f in f_list]
    writing_dict = taxa_acc_dict(dict_list, taxa)
    print "\nWriting accession table."
    write_acc_table(taxa, writing_dict, f_names)
    tf = datetime.now()
    te = tf - tb
    print "\nFinished. Elapsed time: {} (H:M:S)\n\n".format(te)

if __name__ == '__main__':
    main()

