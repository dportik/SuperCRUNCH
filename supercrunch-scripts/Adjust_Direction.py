''''
SuperCRUNCH: Adjust_Direction module

    Adjust_Direction: The purpose of this module is to ensure sequences 
    are all written in the same direction before performing alignments or 
    additional filtering. The input is an unaligned fasta file and 
    the output is an unaligned fasta file with all sequences in the 
    'correct' orientation.
     
    This function processes an unaligned fasta file and adjusts sequence 
    directions by default using the --adjustdirection implementation of 
    mafft. If the optional --accurate flag is included, it will use the 
    --adjustdirectionaccurately option, which is slower but more accurate.
    The number of threads can be specified using the --threads flag.
    The output from mafft is an interleaved fasta with sequences in all 
    lowercase, and sequences that have been reversed are flagged with an 
    '_R_' at the beginning of the record ID. This module takes that file 
    and converts it to a cleaner format. Sequences are written in uppercase, 
    are ungapped (stripping alignment), and the '_R_' is removed from 
    sequence records containing sequences that were reversed. A log file 
    is produced that indicates how many sequences were reversed vs. not, 
    and the record IDs of the sequences requiring reversal are written to 
    a separate log file.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.
    
-------------------------
Compatible with Python 2.7 & 3.7
Python modules required:
	-BioPython
Dependencies:
    -mafft (in path)
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
import argparse
import shutil
import subprocess as sp
from datetime import datetime
from Bio import SeqIO

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Adjust_Direction: The purpose of this module is to ensure sequences 
    are all written in the same direction before performing alignments or 
    additional filtering. The input is an unaligned fasta file and 
    the output is an unaligned fasta file with all sequences in the 
    'correct' orientation.
     
    This function processes an unaligned fasta file and adjusts sequence 
    directions by default using the --adjustdirection implementation of 
    mafft. If the optional --accurate flag is included, it will use the 
    --adjustdirectionaccurately option, which is slower but more accurate.
    The number of threads can be specified using the --threads flag.
    The output from mafft is an interleaved fasta with sequences in all 
    lowercase, and sequences that have been reversed are flagged with an 
    '_R_' at the beginning of the record ID. This module takes that file 
    and converts it to a cleaner format. Sequences are written in uppercase, 
    are ungapped (stripping alignment), and the '_R_' is removed from 
    sequence records containing sequences that were reversed. A log file 
    is produced that indicates how many sequences were reversed vs. not, 
    and the record IDs of the sequences requiring reversal are written to 
    a separate log file.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

    DEPENDENCIES: Python: BioPython; Executables in path: mafft.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i",
                            "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which "
                            "contains the input fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("--threads",
                            default=None,
                            help="OPTIONAL: Specify number of threads to use "
                            "for MAFFT (default is one).")
    
    parser.add_argument("--accurate",
                            action='store_true',
                            help="OPTIONAL: Use --adjustdirectionaccurately MAFFT "
                            "implementation, rather than --adjustdirection.")
    
    return parser.parse_args()

def mafft_adjust(f, accurate, threads, fdir, ldir):
    """
    First generates proper call string to use with sp.call() 
    based on the input file (f) and arguments supplied. Uses
    sp.call() to generate output file, which is then re-written
    as a stripped alignment file using biopython. During the
    re-writing, sequence names are inspected to see if their 
    direction had been changed (denoted by a _R_ at beginning
    of sequence label). Counts of adjusted and non-adjusted
    sequences are generated and a list containing the locus
    name, adjusted count, and non-adjusted count is returned.
    The temporary output file is removed and the stripped 
    alignment is moved to the relevant output directory.
    """
    print("\n\nAdjusting direction of sequences for {}\n\n".format(f))
    tb = datetime.now()
    #find correct filename prefix to use
    prefix = f.split('.')[0]
    if threads is None:
        threads = 1

    #create command line string and use
    if accurate is True:
        call_string = "mafft --thread {0} --adjustdirectionaccurately {1} > {2}_temp.fasta".format(threads, f, prefix)

    else:
        call_string = "mafft --thread {0} --adjustdirection {1} > {2}_temp.fasta".format(threads, f, prefix)
        
    print("{}\n".format(call_string))
    proc = sp.call(call_string, shell=True)
    
    #load adjusted fasta file as indexed dictionary structure
    fasta_dict = SeqIO.index("{0}_temp.fasta".format(prefix), "fasta")
    
    out_fasta, name_log = "{0}_Adjusted.fasta".format(prefix), "{}_Adjusted_Name_Log.txt".format(prefix)
    
    with open(name_log, 'a') as fh_log:
        fh_log.write("Sequences_Flipped\n")
    
    adjusted, fine = int(0), int(0)
    #search for signs a sequence was adjusted, noted by the
    #_R_ at the beginning of the sequence description
    with open(out_fasta, 'a') as fh_out:
        for record in fasta_dict:
            
            #action for sequences that have been reversed ('_R_')
            if fasta_dict[record].description.startswith('_R_'):
                adjusted += 1
                with open(name_log, 'a') as fh_log:
                    fh_log.write("{}\n".format(fasta_dict[record].description))
                    
                new_description = fasta_dict[record].description.replace('_R_', '')
                newseq = fasta_dict[record].seq.upper().ungap("-")
                fh_out.write( ">{}\n{}\n".format(new_description, newseq))
                
            #action for sequences that were not adjusted
            else:
                fine += 1
                newseq = fasta_dict[record].seq.upper().ungap("-")
                fh_out.write( ">{}\n{}\n".format(fasta_dict[record].description, newseq))

    #put locus name and counts into a list to return
    info_list = [prefix, fine, adjusted]

    #cleanup - remove temp file and move other files
    os.remove("{0}_temp.fasta".format(prefix))
    shutil.move(out_fasta, fdir)
    shutil.move(name_log, ldir)

    tf = datetime.now()
    te = tf - tb
    print("Adjusted {} of {} sequences".format(adjusted, fine))
    print("Total time to adjust sequences in {0}: {1} (H:M:S)\n\n".format(f, te))
    
    return info_list

def write_master_log(adjusted_info, rdir):
    """
    Writes main log file which contains the following
    columns:
    Locus, Seqs_Correct_Direction, Seqs_Direction_Adjusted
    The log file is populated with information for all
    loci which is generated using the mafft_adjust() function.
    After writing, log file is moved to the relevant output
    directory.
    """
    log_name = "Sequences_Adjusted.txt"
    
    with open(log_name, 'a') as fh:
        fh.write("Locus\tSeqs_Correct_Direction\tSeqs_Direction_Adjusted\n")

    for a in adjusted_info:
        with open(log_name, 'a') as fh:
            fh.write("{}\t{}\t{}\n".format(a[0], a[1], a[2]))
            
    shutil.move(log_name, rdir)
        

def mafft_adjust_runner(indir, accurate, threads, outdir, fdir, ldir):
    """
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the mafft_adjust function
    for each file found. Moves all output files to the 
    relevant output directories.
    """
    os.chdir(indir)
        
    adjusted_info = []
    
    flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
    
    for f in flist:
        info = mafft_adjust(f, accurate, threads, fdir, ldir)
        adjusted_info.append(info)

    write_master_log(adjusted_info, outdir)

def make_dirs(outdir):
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    os.chdir(outdir)
    #get current path
    curpath = os.getcwd()
    
    #create paths using os.path.join() to avoid any issues
    fdir = os.path.join(curpath, "Adjusted-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)
        
    ldir = os.path.join(curpath, "Log-Files")
    if not os.path.exists(ldir):
        os.mkdir(ldir)
                
    return fdir, ldir

def main():
    tb = datetime.now()
    args = get_args()
    
    fdir, ldir = make_dirs(args.outdir)
    
    print("\n\n--------------------------------------------------------------------------------------")
    print("Beginning sequence direction adjustments using MAFFT...")
    print("--------------------------------------------------------------------------------------")
    
    mafft_adjust_runner(args.indir, args.accurate, args.threads, args.outdir, fdir, ldir)
    tf = datetime.now()
    te = tf - tb
    
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")


if __name__ == '__main__':
    main()
