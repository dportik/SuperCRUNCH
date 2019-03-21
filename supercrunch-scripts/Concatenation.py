'''
SuperCRUNCH: Concatenation module

Usage: python Concatenation.py  -i [full path to directory with all fasta files] (REQUIRED)
                                -f [fasta or phylip] (REQUIRED)
                                -s [dash, N, or ?] (REQUIRED)
                                -o [fasta or phylip] (REQUIRED)
						
    Concatenation: Combine multiple alignments into a single concatenated alignment.
    The input file format can be non-interleaved fasta or phylip, and is selected using the -f flag.
    The input alignment files must be labeled with one of the following 
    extensions to be read: NAME.fasta, NAME.fa, NAME.phylip, or NAME.phy.
    The complete set of taxa is assessed across all alignments. All taxon labels 
    must be unique (no duplicate names) or an error will be thrown. If a taxon is absent
    from an alignment, a missing sequence is generated using the symbol selected with the
    -s flag (can be an N, dash, or ?). The output format must be specified as fasta or phylip
    using the -o flag. 

    This script uses python dictionary structures to read alignments, which allows 
    alignments to be processed and concatenated incredibly fast. It was designed to
    to run on a directory containing thousands of large alignment files. 


    Output files include the following, written to the main directory:

           Data_Partitions.txt - File containing the order of loci concatenated and their
                                 corresponding base pairs within the alignment.
                                 
           Taxa_Loci_Count.log - A simple count of the number of sequences that were available for
                                 each taxon, in other words how many alignments the taxon was present in.
                                  
           Concatenated_Alignment.[fasta or phylip] - The final concatenated alignment in the output
                                                    format specified. Taxa are written in alphabetical order.

-------------------------
For Python 2.7
Dependencies: None
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
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Concatenation: Combine multiple alignments into a single concatenated alignment.
    The input file format can be non-interleaved fasta or phylip, and is selected using the -f flag.
    The input alignment files must be labeled with one of the following 
    extensions to be read: NAME.fasta, NAME.fa, NAME.phylip, or NAME.phy.
    The complete set of taxa is assessed across all alignments. All taxon labels 
    must be unique (no duplicate names) or an error will be thrown. If a taxon is absent
    from an alignment, a missing sequence is generated using the symbol selected with the
    -s flag (can be an N, dash, or ?). The output format must be specified as fasta or phylip
    using the -o flag. 
 
    This script uses python dictionary structures to read alignments, which allows 
    alignments to be processed and concatenated incredibly fast. It was designed to
    to run on a directory containing thousands of large alignment files. 

    Output files include the following, written to the main directory:

           Data_Partitions.txt - File containing the order of loci concatenated and their
                                 corresponding base pairs within the alignment.
           Taxa_Loci_Count.log - A simple count of the number of sequences that were available for
                                 each taxon, in other words how many alignments the taxon was present in. 
           Concatenated_Alignment.[fasta or phylip] - The final concatenated alignment in the output
                                   format specified. Taxa are written in alphabetical order.
    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files.")
    parser.add_argument("-f", "--in_format", required=True, choices=["fasta", "phylip"], help="REQUIRED: The input file format for alignments.")
    parser.add_argument("-s", "--symbol", required=True, choices=["dash", "N", "?"], help="REQUIRED: A base pair symbol used to represent missing data when sequences are not available for a taxon.")
    parser.add_argument("-o", "--out_format", required=True, choices=["fasta", "phylip"], help="REQUIRED: The output file format for the final concatenated alignment.")
    return parser.parse_args()        

def get_taxa_fasta(f):
    '''
    Retrieve the names following the > in all lines
    from a fasta file (f), return as a list.
    '''
    with open(f, 'r') as fh:
        taxa = [l.replace(">",'').replace("\n",'') for l in fh if l.startswith(">")]
    return taxa
            
def get_taxa_phylip(f):
    '''
    Retrieve the names of each sequence
    from a phylip file (f), return as a list.
    '''
    with open(f, 'r') as fh:
        next(fh)
        taxa = [l.split()[0] for l in fh if l.split()]
    return taxa

def collect_taxa(f_list, in_format):
    '''
    Collect taxon names from all alignment files (f_list)
    in specified format (in_format) and return a sorted 
    list of all unique names.
    '''
    taxon_set = set()
    for f in f_list:
        if in_format == "fasta":
            f_taxa = get_taxa_fasta(f)
            taxon_set.update(set(f_taxa))
        elif in_format == "phylip":
            f_taxa = get_taxa_phylip(f)
            taxon_set.update(set(f_taxa))
    taxon_list = sorted(taxon_set)
    return taxon_list

def fasta_dict1(f):
    '''
    ****Not used because fasta formatting should be irrelevant.
    Function to convert non-interleaved fasta file (f) 
    into dictionary structure with taxon as key
    and sequence as value.
    '''
    f_dict = {}
    with open(f, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    for c, line in enumerate(lines):
        if line.startswith(">"):
            f_dict[line.replace(">",'')] = lines[c+1].upper()
    return f_dict

def fasta_dict2(f):
    '''
    Function to convert any fasta file (f) into
    dictionary structure with taxon as key and 
    sequence as value.
    '''
    f_dict = {}
    with open(f, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    for line in lines:
        if line.startswith(">"):
            new_key = line.replace(">",'')
            f_dict[new_key] = ""
        else:
            f_dict[new_key] += line.upper()
    return f_dict
        
def phylip_dict(f):
    '''
    Function to convert phylip file (f) into
    dictionary structure with taxon as key
    and sequence as value.
    '''
    f_dict = {}
    with open(f, 'r') as fh:
        lines = [l.split() for l in fh if l.split()]
    for line in lines[1:]:
        f_dict[line[0].strip()] = line[1].strip().upper()
    return f_dict
    
def collect_dicts(f_list, in_format):
    '''
    Create a list of dictionaries, in which each dictionary is
    created from a list of fasta or phylip files (f_list, in_format),
    with taxon name as the key and sequence as the value. List will 
    be in the same order as the file list (f_list), which was sorted 
    alphabetically.
    '''
    if in_format == "fasta":
        dict_list = [fasta_dict2(f) for f in f_list]
    elif in_format == "phylip":
        dict_list = [phylip_dict(f) for f in f_list]
    return dict_list

def collect_bp_lengths(dict_list):
    '''
    Gather the alignment lengths from the alignment
    dictionaries in the larger list (dict_list). Uses 
    random.choice function which draws an aribitrary 
    key + value pair from the dictionary. 
    List comprehension here gets the length of the 
    value resulting from each random.choice key
    returned for each dictionary, resulting in a 
    list of bp lengths in the same order as the input
    dictionary/input file list.
    '''
    lengths = [len(d[random.choice(d.keys())]) for d in dict_list]
    return lengths

def symbol_dict(symbol):
    '''
    Annoying workaround because argparse hates
    dashes (-) on the cl so we can't use that symbol.
    '''
    sym_dict = {'dash':'-', 'N':'N', '?':'?'}
    sym_val = sym_dict[symbol]
    return sym_val

def create_concat_dict(fdata, taxa, lengths, f_list, sym):
    '''
    Iterate through taxon list (taxa) and create a concatenated
    sequence for each taxon based on alignment dictionaries
    contained in the list (fdata). The lengths list is used to 
    generate an appropriate length of missing data string using
    the symbol selected (sym). The larger dictionary structure
    has taxa as keys and concatenated sequences as their values.
    Counts number of sequences/loci available
    for each taxon and writes to output log file. 
    '''
    concat_dict = {}
    with open("Taxa_Loci_Count.log", 'a') as fh_out:
        fh_out.write("Taxon\tLoci\n")
    for taxon in taxa:
        locus_count = int(0)
        concat_dict[taxon] = ""
        for c, d in enumerate(fdata):
            if taxon in d:
                locus_count += 1
                concat_dict[taxon] += d[taxon]
                #sanity check...
                #print "{} in {}, length {}\n{}\n\n".format(taxon, f_list[c], lengths[c], d[taxon])
            elif taxon not in d:
                missing = int(lengths[c]) * sym
                concat_dict[taxon] += missing
                #print "{} not in {}, length {}\n{}\n\n".format(taxon, f_list[c], lengths[c], missing)
        with open("Taxa_Loci_Count.log", 'a') as fh_out:
            fh_out.write("{}\t{}\n".format(taxon,locus_count))
    return concat_dict

def write_partitions(f_list, lengths):
    '''
    Writes a partitions file that can be used or 
    translated to use with different phylogenetic 
    programs, based on the file name list (f_list)
    and corresponding list of sequence lengths (lengths).
    '''
    bp_count = int(0)
    with open("Data_Partitions.txt", 'a') as fh_out:
        fh_out.write('[data_blocks]\n')
        for c, f in enumerate(f_list):
            begin = bp_count + int(1)
            end = bp_count + int(lengths[c])
            fh_out.write("{0} = {1}-{2};\n".format(f, begin, end))
            bp_count += int(lengths[c])

def write_concatenated(taxa, concat_dict, out_format):
    '''
    Write output file in correct format (out_format) using the 
    concatenated sequence dictionary (concat_dict) and the list of taxa (taxa).
    '''
    print "\tTotal alignment length = {} bp.".format(len(concat_dict[random.choice(concat_dict.keys())]))
    if out_format == "fasta":
        with open("Concatenated_Alignment.fasta", 'a') as fh_out:
            for taxon in taxa:
                fh_out.write(">{0}\n{1}\n".format(taxon, concat_dict[taxon]))
                                
    elif out_format == "phylip":
        with open("Concatenated_Alignment.phy", 'a') as fh_out:
            fh_out.write("{0} {1}\n".format(len(taxa), len(concat_dict[random.choice(concat_dict.keys())])))
            for taxon in taxa:
                fh_out.write("{0} {1}\n".format(taxon, concat_dict[taxon]))

#---------------------------------------------------------------------------------------

def main():
    tb = datetime.now()
    
    args = get_args()
    os.chdir(args.in_dir)
    
    if args.in_format == "fasta":
        f_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])
    elif args.in_format == "phylip":
        f_list = sorted([f for f in os.listdir('.') if f.endswith(".phylip") or f.endswith(".phy")])
    f_list.sort()
    print "\nFound {0} {1} files to concatenate.".format(len(f_list), args.in_format)
    
    taxa = collect_taxa(f_list, args.in_format)
    print "\nFound {} unique taxa across alignment files.".format(len(taxa))
    
    dict_list = collect_dicts(f_list, args.in_format)     
    lengths = collect_bp_lengths(dict_list)
    
    print "\nGathering sequences for all taxa (this could take some time)..."
    sym_val = symbol_dict(args.symbol)
    concat_dict = create_concat_dict(dict_list, taxa, lengths, f_list, sym_val)
    print "\tDone."
    
    print "\nWriting partitions file."
    write_partitions(f_list, lengths)
    
    print "\nWriting concatenated alignment."
    write_concatenated(taxa, concat_dict, args.out_format)
    
    tf = datetime.now()
    te = tf - tb
    print "\nFinished.\nElapsed time: {} (H:M:S)\n".format(te)
    
if __name__ == '__main__':
    main()

