

import os
import shutil
import argparse
import numpy as np
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Seq_Length_Heterogeneity - Calculates several metrics in an effort to quantify
    sequence length variation within alignments. Designed to work for a directory
    (-i) that contains either phylip or fasta alignment files, with the format
    specified by the -f flag. The user must select the sequence type
    to calculate metrics from (-s), with choices 'bases' or 'alns'. For 'bases', 
    only the nucleotides are counted for a sequence (all gaps are ignored), and the
    length is obtained from the sum of the number of bases present. For
    'alns', the position of the first nucleotide and last nucleotide in 
    the sequence alignment are used to infer the aligned sequence length. 
    From the lengths obtained from all the sequences in an alignment, the 
    following metrics are calculated: mean, minimum, maximum, standard deviation, 
    and coefficient of variation. The results for all alignments are written to 
    a tab-delimited output file in the output directory specified (-o). 
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--in_dir",
                            required=True,
                            help="REQUIRED: The full path to a directory which "
                            "contains the input alignment files (fasta or phylip "
                            "format).")
    
    parser.add_argument("-o", "--out_dir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-f", "--format",
                            required=True,
                            choices=["fasta", "phylip"],
                            help="REQUIRED: The file format of the INPUT alignments.")
    
    parser.add_argument("-s", "--seqtype",
                            required=True,
                            choices=["bases", "alns"],
                            help="REQUIRED: Use de-aligned sequence lengths (bases) "
                            "or aligned sequence lengths (alns) for calculations.")
    
    return parser.parse_args()        

def fasta_dict(f):
    """
    Function to convert any fasta file (f) into
    dictionary structure with taxon as key and 
    sequence as value.
    """
    fdict = {}
    with open(f, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
        
    for line in lines:
        if line.startswith(">"):
            new_key = line.replace(">",'')
            fdict[new_key] = ""
            
        else:
            fdict[new_key] += line.upper()
            
    return fdict
        
def phylip_dict(f):
    """
    Function to convert phylip file (f) into
    dictionary structure with taxon as key
    and sequence as value.
    """
    fdict = {}
    
    with open(f, 'r') as fh:
        lines = [l.split() for l in fh if l.split()]
        
    for line in lines[1:]:
        fdict[line[0].strip()] = line[1].strip().upper()
        
    return fdict

def base_calc(seq):
    """
    Function to count the number of bases in
    any sequence, as long as the character
    is not a - or ?. Returns count, which is
    equivalent to the length of the seq, which
    is an integer.
    """
    bcount = sum(1 for b in seq if b != "-" and b != "?")
    
    return bcount

def aln_calc(seq):
    """
    Function to calculate the length of
    an ALIGNED sequence. Finds the first and
    last nucleotide in the aligned sequence 
    and calculates length from the difference
    in these positions. Returns length, which
    is an integer.
    """
    indices = []
    bases = ["A", "C", "T", "G",
                 "R", "Y", "S", "W",
                 "K", "M", "B", "D",
                 "H", "V"]
    for b in bases:
        try:
            i = seq.index(b)
            indices.append(i)
            j = seq.rindex(b)
            indices.append(j)
            #print "\t\t{}: {}, {}".format(b,i,j)
            
        except:
            pass
        
    indices.sort()
    length = indices[-1] - indices[0]
    #print "{} - {} = {} length".format(indices[-1], indices[0], length)
    
    return length

def get_lengths(fdict, seqtype):
    """
    Function to iterate over sequences in a sequence
    dictionary and create a list of all the lengths.
    Returns the list of lengths, which are all integers.
    """
    if seqtype == "bases":
        lengths = [base_calc(fdict[rec]) for rec in fdict]
        
    elif seqtype == "alns":
        lengths = [aln_calc(fdict[rec]) for rec in fdict]
        
    return lengths

def get_stats(l):
    """
    Function to calculate basic stats for the 
    list of lengths, which are all integers.
    """
    avg = np.mean(l)
    std = np.std(l)
    cov = std / avg
    mi = np.amin(l)
    ma = np.amax(l)
    stats = [np.round(avg, 1), mi, ma, np.round(std, 1), np.round(cov, 3)]
    
    return stats

def run_calcs(flist, form, seqtype):
    """
    Main function. Converts input file into a
    sequence dictionary. Counts the lengths of
    all sequences in the sequence dictionary based
    on argument seqtype (e.g. aligned or unaligned
    sequences). Appends information to results list,
    and returns results list.
    """
    results = []
    
    for f in flist:
        fdict = fasta_dict(f) if form == "fasta" else phylip_dict(f)
        lengths = get_lengths(fdict, seqtype)
        stats = get_stats(lengths)
        results.append([f, stats])
        print("\n{0}:\n\t{1}".format(f, stats))
        
    return results

def write_output(results, out_dir, seqtype):
    """
    Function to write output file of information
    which comes from the results list.
    """
    os.chdir(out_dir)
    
    if seqtype == "bases":
        outname = "Sequence_Length_Heterogeneity_Results_BasePairLengths.txt"
        
    elif seqtype == "alns":
        outname = "Sequence_Length_Heterogeneity_Results_AlignedLengths.txt"
        
    with open(outname, 'a') as fh:
        fh.write("Alignment_File\tMeanLength\tMinLength\tMaxLength\tSD\tCoV\n")
        for r in results:
            fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
                         .format(r[0], r[1][0], r[1][1], r[1][2], r[1][3], r[1][4]))
    
def main():
    tb = datetime.now()
    
    args = get_args()
    os.chdir(args.in_dir)
    
    if args.format == "fasta":
        flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
            
    elif args.format == "phylip":
        flist = sorted([f for f in os.listdir('.') if f.endswith((".phylip", ".phy"))])
        
    if not flist:
        print("\nNo files matching the format specified ({}) were found.".format(args.format))
        print("Please check if input file format (-f) was specified correctly.")
        print("If so, please ensure the files contain the proper file extensions: ")
        print("\tfasta - .fasta or .fa\n\tphylip - .phy or .phylip\n")
        
    else:
        results = run_calcs(flist, args.format, args.seqtype)
        if not results:
            print("\nNo results generated, something went wrong.")
            print("Please check file formats are correct.\n")
        else:    
            write_output(results, args.out_dir, args.seqtype)
    
            tf = datetime.now()
            te = tf - tb
            print("\n\n--------------------------------------------------------------------------------------")
            print("\nFinished. Elapsed time: {0} (H:M:S)\n".format(te))
            print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()

