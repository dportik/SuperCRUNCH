'''
SuperCRUNCH: Concatenation module
						
    Concatenation: Combine multiple alignments into a single concatenated 
    alignment. The input file format can be fasta or phylip, 
    and is selected using the --informat flag. The input alignment files must be 
    labeled with one of the following extensions to be read: NAME.fasta, 
    NAME.fa, NAME.phylip, or NAME.phy. The complete set of taxa is assessed 
    across all alignments. All taxon labels must be unique within each alignment
    (no duplicate names) or an error will be thrown. If a taxon is absent 
    from an alignment, a missing sequence is generated using the symbol 
    selected with the -s flag (can be an N, dash, or ?). The output format 
    must be specified as fasta or phylip using the --outformat flag. 

    Output files are written to the specified output directory (-o).

-------------------------
Compatible with Python 2.7 & 3.7
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
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Concatenation: Combine multiple alignments into a single concatenated 
    alignment. The input file format can be fasta or phylip, 
    and is selected using the --informat flag. The input alignment files must be 
    labeled with one of the following extensions to be read: NAME.fasta, 
    NAME.fa, NAME.phylip, or NAME.phy. The complete set of taxa is assessed 
    across all alignments. All taxon labels must be unique within each alignment
    (no duplicate names) or an error will be thrown. If a taxon is absent 
    from an alignment, a missing sequence is generated using the symbol 
    selected with the -s flag (can be an N, dash, or ?). The output format 
    must be specified as fasta or phylip using the --outformat flag. 

    Output files are written to the specified output directory (-o).

    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input alignment files (fasta or phylip).")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to "
                            "write output files.")
    
    parser.add_argument("--informat",
                            required=True,
                            choices=["fasta", "phylip"],
                            help="REQUIRED: The file format of the INPUT alignments.")

    parser.add_argument("--outformat",
                            required=True,
                            choices=["fasta", "phylip", "nexus"],
                            help="REQUIRED: The file format for the OUTPUT concatenated "
                            "alignment.")
    
    parser.add_argument("-s", "--symbol",
                            required=True,
                            choices=["dash", "N", "?"],
                            help="REQUIRED: A base pair symbol used to represent missing "
                            "data when sequences are not available for a taxon.")
    
    return parser.parse_args()        

def get_taxa_fasta(f):
    """
    Retrieve the names following the > in all lines
    from a fasta file (f), return as a list.
    """
    with open(f, 'r') as fh:
        taxa = set([l.replace(">",'').replace("\n",'') for l in fh
                    if l.startswith(">")])
            
    return taxa
            
def get_taxa_phylip(f):
    """
    Retrieve the names of each sequence
    from a phylip file (f), return as a set.
    """
    with open(f, 'r') as fh:
        #must skip first line which has sequence and bp info
        next(fh)
        #for remaining lines use list comprehension to
        #split lines to obtain taxon names and put in list
        taxa = set([l.split()[0] for l in fh if l.split()])
        
    return taxa

def collect_taxa(flist, informat):
    """
    Collect taxon names from all alignment files (flist)
    in specified format (in_format). Each file is run 
    through a get_taxa function, and the names 
    returned as a temporary set are added to the larger 
    set. The larger set is converted to a sorted list
    and returned.
    """
    taxon_set = set()
    
    for f in flist:
        if informat == "fasta":
            ftaxa = get_taxa_fasta(f)
            taxon_set.update(ftaxa)
            
        elif informat == "phylip":
            ftaxa = get_taxa_phylip(f)
            taxon_set.update(ftaxa)
            
    taxon_list = sorted(taxon_set)
    
    return taxon_list

def fasta_dict(f):
    """
    Function to convert any fasta file (f) into
    a dictionary structure with taxon as key and 
    sequence as value. Avoids biopython to read.
    """
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
    """
    Function to convert phylip file (f) into
    dictionary structure with taxon as key
    and sequence as value. Avoids biopython 
    to read files.
    """
    f_dict = {}
    
    with open(f, 'r') as fh:
        lines = [l.split() for l in fh if l.split()]
        
    for line in lines[1:]:
        f_dict[line[0].strip()] = line[1].strip().upper()
        
    return f_dict
    
def collect_dicts(flist, informat):
    """
    Create a list of dictionaries, in which each dictionary is
    created from a list of fasta or phylip files (flist, in_format),
    with taxon name as the key and sequence as the value. List will 
    be in the same order as the file list (flist), which was sorted 
    alphabetically.
    """
    if informat == "fasta":
        dict_list = [fasta_dict(f) for f in flist]
        
    elif informat == "phylip":
        dict_list = [phylip_dict(f) for f in flist]
        
    return dict_list

def collect_bp_lengths(dict_list):
    """
    Gather the alignment lengths from the alignment
    dictionaries in the larger list (dict_list). Uses 
    random.choice function which draws an aribitrary 
    key from the dictionary, and obtains the length
    of the value (aligned sequence) based on the key.
    List comprehension here gets the length of the 
    value resulting from each random.choice key
    returned for each dictionary, resulting in a 
    list of bp lengths in the same order as the input
    dictionary/input file list.
    """
    lengths = [len(d[random.choice(list(d.keys()))]) for d in dict_list]
    
    return lengths

def symbol_dict(symbol):
    """
    Workaround because argparse hates
    dashes (-) on the cl so we can't use that symbol.
    Instead, we lookup a value in a dictionary and
    return that.
    """
    sym_dict = {'dash':'-', 'N':'N', '?':'?'}
    sym_val = sym_dict[symbol]
    
    return sym_val

def get_seq_numbers(taxloci):
    with open(taxloci, 'r') as fh:
        next(fh)
        seqs = [int(l.split()[1]) for l in fh if l.strip()]
    seq_count = sum(seqs)

    return seq_count
    
def create_concat_dict(fdata, taxa, lengths, flist, sym):
    """
    Iterate through taxon list (taxa) and create a concatenated
    sequence for each taxon based on alignment dictionaries
    contained in the list (fdata). The lengths list is used to 
    generate an appropriate length of missing data string using
    the symbol selected (sym). The larger dictionary structure
    has taxa as keys and concatenated sequences as their values.
    Counts number of sequences/loci available
    for each taxon and writes to output log file. 
    """
    concat_dict = {}
    taxloci = "Taxa_Loci_Count.log"
    
    with open(taxloci, 'a') as fh:
        fh.write("Taxon\tLoci\n")
        
    for taxon in taxa:
        locus_count = int(0)
        concat_dict[taxon] = ""
        #note enumerate was only necessary for testing
        #with the hashed out print statements
        for c, d in enumerate(fdata):
            if taxon in d:
                locus_count += 1
                concat_dict[taxon] += d[taxon]
                #sanity check...
                #print("{} in {}, length {}\n{}\n\n".format(taxon, flist[c], lengths[c], d[taxon]))
                
            elif taxon not in d:
                missing = int(lengths[c]) * sym
                concat_dict[taxon] += missing
                #print("{} not in {}, length {}\n{}\n\n".format(taxon, flist[c], lengths[c], missing))
                
        with open(taxloci, 'a') as fh:
            fh.write("{}\t{}\n".format(taxon, locus_count))

    seq_count = get_seq_numbers(taxloci)
            
    return concat_dict, seq_count

def write_partitions(flist, lengths):
    """
    Writes a partitions file that can be used or 
    translated to use with different phylogenetic 
    programs, based on the file name list (flist)
    and corresponding list of sequence lengths (lengths).
    """
    bp_count = int(0)
    
    with open("Data_Partitions.txt", 'a') as fh_out:
        fh_out.write('[data_blocks]\n')
        
        for c, f in enumerate(flist):
            begin = bp_count + int(1)
            end = bp_count + int(lengths[c])
            fh_out.write("{0} = {1}-{2};\n".format(f, begin, end))
            bp_count += int(lengths[c])

def write_concatenated(taxa, concat_dict, outformat, symbol):
    """
    Write output file in correct format (out_format) using the 
    concatenated sequence dictionary (concat_dict) and the list of taxa (taxa).
    """    
    if outformat == "fasta":
        
        with open("Concatenated_Alignment.fasta", 'a') as fh:
            for taxon in taxa:
                fh.write(">{0}\n{1}\n".format(taxon, concat_dict[taxon]))
                                
    elif outformat == "phylip":
        
        with open("Concatenated_Alignment.phy", 'a') as fh:
            fh.write("{0} {1}\n".format(len(taxa), len(concat_dict[random.choice(list(concat_dict.keys()))])))
            for taxon in taxa:
                fh.write("{0} {1}\n".format(taxon, concat_dict[taxon]))
                
    elif outformat == "nexus":
        with open("Concatenated_Alignment.nex", 'a') as fh:
            fh.write('''
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=DNA  MISSING={2} GAP=-;
MATRIX
'''.format(len(taxa), len(concat_dict[random.choice(list(concat_dict.keys()))]), symbol))
            for taxon in taxa:
                fh.write("\n{0} {1}".format(taxon, concat_dict[taxon]))
            fh.write("\n;\nEnd;")

def main():
    tb = datetime.now()
    
    args = get_args()
    os.chdir(args.indir)
    
    if args.informat == "fasta":
        flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
        
    elif args.informat == "phylip":
        flist = sorted([f for f in os.listdir('.') if f.endswith((".phylip", ".phy"))])
        
    print("\nFound {:,} {} files to concatenate.".format(len(flist), args.informat))
    
    taxa = collect_taxa(flist, args.informat)
    print("\nFound {:,} unique taxa across alignment files.".format(len(taxa)))
    
    dict_list = collect_dicts(flist, args.informat)     
    lengths = collect_bp_lengths(dict_list)

    os.chdir(args.outdir)
    
    print("\nGathering sequences for all taxa (this could take some time)...")
    sym_val = symbol_dict(args.symbol)
    concat_dict, seq_count = create_concat_dict(dict_list, taxa, lengths, flist, sym_val)
    print("\tDone.")
    
    print("\nWriting partitions file.")
    write_partitions(flist, lengths)
    
    print("\nWriting concatenated alignment.")
    print("\tTotal alignment length = {:,} bp."
              .format(len(concat_dict[random.choice(list(concat_dict.keys()))])))
    print("\tTotal number of sequences included = {:,}.".format(seq_count))

    write_concatenated(taxa, concat_dict, args.outformat, sym_val)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
if __name__ == '__main__':
    main()

