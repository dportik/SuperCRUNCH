'''
SuperCRUNCH: Make_Acc_Table module

    Make_Acc_Table: A tool for extracting taxon and accession number 
    information from a directory of fasta files containing records 
    with the original description lines.
    For example:
    >FJ745388.1 Abelia chinensis voucher N. Pyck 1989-2220 internal ...

    A comprehensive list of taxa is collected across all fasta files 
    in the directory. The taxon names to search are binomial by default 
    (only genus + species), but if the optional -s flag is used then 
    the subspecies names contained in the supplied text file will also 
    be searched for and included (if present). The table is written with 
    species in the first column (sorted alphabetically) and the fasta 
    files composing the remaining columns (labeled by fasta name), with 
    rows filled with accession numbers or dashes depending on the presence 
    or absence of the taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files. 

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

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
import os
import shutil
import argparse
import random
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Make_Acc_Table: A tool for extracting taxon and accession number 
    information from a directory of fasta files containing records 
    with the original description lines.
    For example:
    >FJ745388.1 Abelia chinensis voucher N. Pyck 1989-2220 internal ...

    A comprehensive list of taxa is collected across all fasta files 
    in the directory. The taxon names to search are binomial by default 
    (only genus + species), but if the optional -s flag is used then 
    the subspecies names contained in the supplied text file will also 
    be searched for and included (if present). The table is written with 
    species in the first column (sorted alphabetically) and the fasta 
    files composing the remaining columns (labeled by fasta name), with 
    rows filled with accession numbers or dashes depending on the presence 
    or absence of the taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files. 

    To be recognized, a fasta file must be labeled as "NAME.fasta" or 
    "NAME.fa". The NAME portion should not contain any periods or spaces, 
    but can contain underscores. 
    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-s", "--subspecies",
                            required=False,
                            default=None,
                            help="OPTIONAL: The full path to a text file containing all "
                            "subspecies names to cross-reference in the fasta file.")
    
    return parser.parse_args()        

def parse_taxa(f):
    """
    Retrieve subspecies (trinomial) names ONLY from user supplied taxon 
    names file (f). Will return sorted list of subspecies, with all 
    names in uppercase.
    """
    print("\nParsing taxon information from {}.".format(f))
    
    with open(f, 'r') as fh:
        subspecies_set = set([line.upper().strip().replace(" ","_") for line in fh
                                  if len(line.split()) == int(3)])
        
    subspecies = sorted(subspecies_set)
    print("\tFound {} subspecies names to include.\n".format(len(subspecies)))
    
    return subspecies

def get_taxon(line):
    """
    Retrieve the taxon name from '>' line (line) in a fasta file.
    Will fetch the species (binomial) and subspecies (trinomial)
    names and return both names. These names will be in uppercase.
    """
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if line.split() >= int(3)][1:3]
    taxon_sp = "_".join(parts1)
    
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if line.split() >= int(4)][1:4]
    taxon_ssp = "_".join(parts2)
    
    return taxon_sp, taxon_ssp

def get_taxa_fasta(f, subspecies):
    """
    Retrieve the names following the > in all lines
    from a fasta file (f) based on subspecies option
    (subspecies = True or False), add names to the set,
    then return the set as a sorted list.
    """
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
    
def collect_taxa(flist, subspecies):
    """
    Collect taxon names from all alignment files (flist)
    based on subspecies option. Create a set of names 
    from each file, add to a larger set, then return 
    the larger set as a sorted list.
    """
    taxon_set = set()
    for f in flist:
        ftaxa = get_taxa_fasta(f, subspecies)
        taxon_set.update(set(ftaxa))
    taxon_list = sorted(taxon_set)
    
    return taxon_list

def f_taxon_dict(f, subspecies):
    """
    Function to convert fasta file (f) into
    dictionary structure with taxon label as key
    and the accession number as the value, based 
    on the subspecies option.
    """
    fdict = {}
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
                    fdict[taxon_ssp] = acc
                    
                else:
                    fdict[taxon_sp] = acc
                    
            else:
                fdict[taxon_sp] = acc
                
    return fdict

def taxa_acc_dict(dict_list, taxa):
    """
    Create a dictionary where keys are taxa from
    taxa list (taxa) and values are sequences from the
    alignments (or are supplied as missing data values). 
    All fasta files were previously converted to dictionaries 
    with taxa as keys and seqs as values, then put into 
    a larger list (dict_list) which is in alphabetical
    order by alignment file name (that way we know which
    alignment a seq is coming from across the list).
    """
    writing_dict = {}
    for taxon in taxa:
        writing_dict[taxon] = ""
        for d in dict_list:
            if taxon in d:
                writing_dict[taxon] += "{}\t".format(d[taxon])
            elif taxon not in d:
                writing_dict[taxon] += "{}\t".format("-")
                
    return writing_dict

def write_acc_table(taxa, writing_dict, fnames):
    """
    Writes an accession table for all entries of the taxon list (taxa).
    Uses the dictionary structure (writing_dict) in which keys are taxa
    and values are the previously concatenated sequences. Labels columns
    appropriately with alignment file names (f_names).
    """
    with open("GenBank_Accession_Table.txt", 'a') as fh:
        fh.write('Taxon\t')
        for f in fnames:
            fh.write('{}\t'.format(f))
        fh.write('\n')
        for taxon in taxa:
            fh.write("{0}\t{1}\n".format(taxon.capitalize().replace("_"," "),
                                                 writing_dict[taxon].strip('\t')))

def main():
    tb = datetime.now()
    args = get_args()
    subspecies = None
    if args.subspecies is not None:
        subspecies = parse_taxa(args.subspecies)
        
    os.chdir(args.indir)
    flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
    print("\nFound {0} files to use for accession table.".format(len(flist)))
    fnames = [f.split(".")[0] for f in flist]
    
    taxa = collect_taxa(flist, subspecies)
    print("\nFound {} unique taxa across fasta files.".format(len(taxa)))
    print("\nGathering accession numbers for all taxa.")
    dict_list = [f_taxon_dict(f, subspecies) for f in flist]
    writing_dict = taxa_acc_dict(dict_list, taxa)

    os.chdir(args.outdir)
    print("\nWriting accession table.")
    write_acc_table(taxa, writing_dict, fnames)
    
    tf = datetime.now()
    te = tf - tb
    print("\nFinished. Elapsed time: {} (H:M:S)\n\n".format(te))

if __name__ == '__main__':
    main()

