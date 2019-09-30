'''
SuperCRUNCH: Filter_Fasta_by_Min_Seqs module

    Filter_Fasta_by_Min_Seqs: The purpose of this module is to examine 
    fasta files for the number of sequences contained within. If taxa 
    are represented by a single sequence, this is equivalent to the 
    number of taxa. Running the script with only the -i and --min_seqs 
    flags will simply print the number of fasta files passing and failing the 
    minimum sequence filter. The minimum number of sequences is an integer 
    specified using the --min_seqs flag. Running this module with the 
    -i, --min_seqs, and the optional -o flags will copy the fasta files 
    that pass the minimum sequence filter to the output directory 
    specified with the -o flag.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores.

-------------------------
Compatible with Python 2.7 & 3.7
Dependencies: 
    None
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
February 2019
Distributed under the 
GNU General Public Lincense
'''
import os
import argparse
import shutil
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Filter_Fasta_by_Min_Seqs: The purpose of this module is to examine 
    fasta files for the number of sequences contained within. If taxa 
    are represented by a single sequence, this is equivalent to the 
    number of taxa. Running the script with only the -i and --min_seqs 
    flags will simply print the number of fasta files passing and failing the 
    minimum sequence filter. The minimum number of sequences is an integer 
    specified using the --min_seqs flag. Running this module with the 
    -i, --min_seqs, and the optional -o flags will copy the fasta files 
    that pass the minimum sequence filter to the output directory 
    specified with the -o flag.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

    DEPENDENCIES: None.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory of fasta files.")
    
    parser.add_argument("--min_seqs",
                            type=int,
                            required=True,
                            help="REQUIRED: The minimum number of taxa/sequences required.")
    
    parser.add_argument("-o", "--outdir",
                            required=False,
                            help="OPTIONAL: The full path to an existing directory to copy "
                            "files passing the filter.")
    
    return parser.parse_args()

def evaluate_recs(f, min_taxa):
    """
    Function to count number of sequence records and evaluate
    against the minimum number requirement supplied by user.
    Returns as True or False.
    """
    passed = False
    
    with open(f, 'r') as fh:
        count = sum(1 for line in fh if line.startswith(">"))
        
    if count >= min_taxa:
        passed = True
        
    return passed   

def filter_fasta(indir, outdir, min_seqs):
    """
    Move to input directory, identify fasta files, then
    evaluate each file using the evaluate_recs() function. 
    If the file passes, copy it to the output directory
    supplied by the user.
    """
    os.chdir(indir)
    
    pass_count, fail_count = int(0), int(0)
    
    flist = [f for f in os.listdir('.') if f.endswith(('.fa', '.fasta'))]
    print("\n\nFound {0:,} total fasta files.\nProcessing...".format(len(flist)))
    
    for f in flist:
        val = evaluate_recs(f, min_seqs)
        
        if val is True:
            shutil.copy(f, outdir)
            pass_count += 1
            
        elif val is False:
            fail_count += 1
            
    print("\n\n{0:,} fasta files passed the minimum sequence filter (>={2} seqs)."
              "\n{1:,} fasta files failed the minimum sequence filter (>={2} seqs).\n\n"
              .format(pass_count, fail_count, min_seqs))
    
    print("Fasta files with >= {0} records have been copied to: {1}\n\n"
              .format(min_seqs, outdir))

def count_fasta(indir, min_seqs):
    """
    Move to input directory, identify fasta files, then
    evaluate each file using the evaluate_recs() function. 
    Rather than copy files, simply print the results.
    """
    os.chdir(indir)
    
    pass_count, fail_count = int(0), int(0)
    
    flist = [f for f in os.listdir('.') if f.endswith(('.fa', '.fasta'))]
    print("\n\nFound {0:,} total fasta files.\nProcessing...".format(len(flist)))
    
    for f in flist:
        val = evaluate_recs(f, min_seqs)
        
        if val is True:
            pass_count += 1
            
        elif val is False:
            fail_count += 1
            
    print("\n\n{0:,} fasta files passed the minimum sequence filter (>={2} seqs)."
              "\n{1:,} fasta files failed the minimum sequence filter (>={2} seqs).\n\n"
              .format(pass_count, fail_count, min_seqs))
    
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
    fdir = os.path.join(curpath, "Filtered-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)

    return fdir
        
def main():
    args = get_args()
    tb = datetime.now()
    
    if args.outdir:
        fdir = make_dirs(args.outdir)
        filter_fasta(args.indir, fdir, args.min_seqs)
        
    else:
        count_fasta(args.indir, args.min_seqs)
        
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Elapsed time : {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
if __name__ == '__main__':
    main()
