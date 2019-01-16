''''
SuperCRUNCH: Adjust_Direction module

Usage: python Adjust_Direction.py -i [directory with all fasta files] (REQUIRED)
                                  --accurate (OPTIONAL)

    Adjust_Direction: The purpose of this script is to check sequences to ensure their proper direction
    before performing alignments or additional filtering (for example based on translation). 
    The input is an unaligned fasta file and the output is an unaligned fasta file with 
    all sequences in the 'correct' direction.
     
    This function processes an unaligned fasta file and adjusts sequence directions by default 
    using the --adjustdirection implementation of mafft. If the optional --accurate flag is included,
    it will use the --adjustdirectionaccurately option, which is slower but more accurate.
    The output from mafft is an interleaved fasta with sequences in all lowercase, and sequences 
    that have been reversed are flagged with an '_R_' at the beginning of the record ID. This 
    script takes that file and converts it to a cleaner format. Sequences are written in
    uppercase, are ungapped (stripping alignment), and the '_R_' is removed from
    sequence records containing sequences that were reversed. A log file is produced
    that indicates how many sequences were reversed vs. not, and the record IDs of the
    sequences requiring reversal are written to a separate log file.

    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.
    
    A few output files are produced per fasta and moved to the output directory called 
    Output_mafft_Adjust_Fasta_Files:
    
    	[NAME]_Direction.fasta - The UNALIGNED fasta file from mafft adjustment,
                                 in which all sequences are now correctly oriented.
    	
    	[NAME]_Adjusted_Name_Log.txt - Contains the full names of the sequences 
                                       that were reversed in this particular fasta file.
    
    Another single output file is created:
    
        Log_Sequences_Adjusted.txt - Contains the names of all fasta files and the number of 
                                     sequences that were found correctly oriented or had to 
                                     be reversed.


-------------------------
Written for Python 2.7
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
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Adjust_Direction: The purpose of this script is to check sequences to ensure their proper direction
    before performing alignments or additional filtering (for example based on translation). 
    The input is an unaligned fasta file and the output is an unaligned fasta file with 
    all sequences in the 'correct' direction.
     
    This function processes an unaligned fasta file and adjusts sequence directions by default 
    using the '--adjustdirection' implementation of mafft. If the optional --accurate flag is included,
    it will use the '--adjustdirectionaccurately' option, which is slower but more accurate.
    The output from mafft is an interleaved fasta with sequences in all lowercase, and sequences 
    that have been reversed are flagged with an '_R_' at the beginning of the record ID. This 
    script takes that file and converts it to a cleaner format. Sequences are written in
    uppercase, are ungapped (stripping alignment), and the '_R_' is removed from
    sequence records containing sequences that were reversed. A log file is produced
    that indicates how many sequences were reversed vs. not, and the record IDs of the
    sequences requiring reversal are written to a separate log file.
    
    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

    
    A few output files are produced per fasta and moved to the output directory called 
    Output_mafft_Adjust_Fasta_Files:
    
    	[NAME]_Direction.fasta - The UNALIGNED fasta file from mafft adjustment, in which all sequences are now correctly oriented.
    	    	
    	[NAME]_Adjusted_Name_Log.txt - Contains the full names of the sequences that were reversed in this particular fasta file.
    
    Another single output file is created:
    
        	Log_Sequences_Adjusted.txt - Contains the names of all fasta files and the number of 
    								sequences that were found correctly oriented or had to 
    								be reversed.
    DEPENDENCIES: Python: BioPython; Executables in path: mafft.
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files. Follow labeling format: NAME.fasta")
    parser.add_argument("--accurate", action='store_true', help="OPTIONAL: Use --adjustdirectionaccurately MAFFT implementation, rather than --adjustdirection.")
    return parser.parse_args()

def directory_mafft_adjust(in_dir, accurate):
    '''
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the mafft_adjust function
    for each file found. Moves all output files to the output
    directory: /Output_mafft_Adjust_Fasta_Files
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "Beginning sequence direction adjustments using MAFFT..."
    print "--------------------------------------------------------------------------------------"
    os.chdir(in_dir)

    out_dir = 'Output_mafft_Adjust_Fasta_Files'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    log_name = "Log_Sequences_Adjusted.txt"
    with open(log_name, 'a') as fh_log:
        fh_log.write("Locus\tSeqs_Correct_Direction\tSeqs_Direction_Adjusted\n")

    f_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])
    for f in f_list:
        summary = mafft_adjust(f, accurate)
        with open(log_name, 'a') as fh_log:
            fh_log.write("{}\t{}\t{}\n".format(summary[0],summary[1],summary[2]))
        for fout in os.listdir('.'):
            if fout.endswith('Adjusted.fasta') or fout.endswith("_Adjusted_Name_Log.txt"):
                shutil.move(fout, out_dir)               
    shutil.move(log_name, out_dir)
    print "\n\n--------------------------------------------------------------------------------------"
    print "Finished sequence direction adjustments."
    print "--------------------------------------------------------------------------------------\n\n"

def mafft_adjust(f, accurate):
    print "\n\nAdjusting direction of sequences for {}\n\n".format(f)
    tb = datetime.now()
    #find correct filename prefix to use
    prefix = f.split('.')[0]
    #create command line string and use
    if accurate is True:
        call_string = "mafft --adjustdirectionaccurately {0} > {1}_temp.fasta".format(f, prefix)
    else:
        call_string = "mafft --adjustdirection {0} > {1}_temp.fasta".format(f, prefix)
    print call_string, '\n'
    proc = sp.call(call_string, shell=True)
    #load adjusted fasta file as indexed dictionary structure
    fasta_dict = SeqIO.index("{0}_temp.fasta".format(prefix), "fasta")
    
    out_fasta = "{0}_Adjusted.fasta".format(prefix)
    name_log = "{}_Adjusted_Name_Log.txt".format(prefix)
    with open(name_log, 'a') as fh_log:
        fh_log.write("Sequences_Flipped\n")
    
    adjusted = int(0)
    fine = int(0)
    #search for signs a sequence was adjusted, noted by the _R_ at the beginning of the seq description
    #record adjusted or not and ungap all seqs to strip alignments
    with open(out_fasta, 'a') as fh_out:
        for record in fasta_dict:
            if fasta_dict[record].description.startswith('_R_'):
                adjusted+=1
                with open(name_log, 'a') as fh_log:
                    fh_log.write("{}\n".format(fasta_dict[record].description))
                new_description = fasta_dict[record].description.strip('_R_')
                newseq = fasta_dict[record].seq.upper().ungap("-")
                fh_out.write( ">{}\n{}\n".format(new_description, newseq))
            else:
                fine+=1
                newseq = fasta_dict[record].seq.upper().ungap("-")
                fh_out.write( ">{}\n{}\n".format(fasta_dict[record].description, newseq))
            
    os.remove("{0}_temp.fasta".format(prefix))

    tf = datetime.now()
    te = tf - tb
    print "Total time to adjust sequences in {0}: {1} (H:M:S)\n\n".format(f,te)
    
    info_list = [prefix,fine,adjusted]
    return info_list
    
#-----------------------------------------------------------------------------------------

def main():
    tb = datetime.now()
    args = get_args()
    directory_mafft_adjust(args.in_dir, args.accurate)
    tf = datetime.now()
    te = tf - tb
    print "\n\nTotal time to adjust sequences across all fasta files: {0} (H:M:S)\n\n".format(te)

if __name__ == '__main__':
    main()
