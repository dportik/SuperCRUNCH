''''
SuperCRUNCH: Fasta_Convert module

    Fasta_Convert: Converts a directory of aligned fasta files into 
    both phylip and nexus formats. Moves new files into corresponding 
    directories in the specified output directory. Note this should 
    only be used after records have been renamed with species or 
    accession names, as the original longer description
    lines containing whitespace will cause problems.
    
    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

    Output_Phylip_Files:
           [NAME].phy - The phylip format alignment.
                                           
    Output_Nexus_Files:
           [NAME].nex - The nexus format alignment
           
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

import argparse
import os
import shutil
from datetime import datetime
from Bio import SeqIO

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Fasta_Convert: Converts a directory of aligned fasta files into 
    both phylip and nexus formats. Moves new files into corresponding 
    directories in the specified output directory. Note this should 
    only be used after records have been renamed with species or 
    accession names, as the original longer description
    lines containing whitespace will cause problems.
    
    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME. 

    DEPENDENCIES: Python: BioPython.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the fasta format alignments to convert.")
    
    parser.add_argument("-o", "--outdir",
                            required=False,
                            help="OPTIONAL: The full path to an existing directory to write "
                            "output files.")

    return parser.parse_args()


def read_fasta(f):
    """
    Quickly obtain the contents of a fasta file in
    biopython list format. Obtain the number of records
    and length of the alignment from this list. Return
    both lists.
    """
    contents = list(SeqIO.parse(f, 'fasta'))
    aln_info = [len(contents), len(contents[0].seq)]            
    return contents, aln_info

def convert_to_phy(contents, aln_info, prefix, pdir):
    """
    Function to write sequence records in phylip format.
    """
    print("\tConverting {}.fasta into phylip format".format(prefix))
    out = "{}.phy".format(prefix)
    
    with open(out, 'a') as fh_out:
        fh_out.write("{0} {1}".format(aln_info[0], aln_info[1]))
        for rec in contents:
            fh_out.write("\n{} {}".format(rec.description, str(rec.seq)))
            
    shutil.move(out, pdir)
    
def convert_to_nex(contents, aln_info, prefix, ndir):
    """
    Function to write sequence records in nexus format.
    """
    print("\tConverting {}.fasta into nexus format\n".format(prefix))
    out = "{}.nex".format(prefix)
    with open(out, 'a') as fh:
        fh.write('''
#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=DNA  MISSING=N GAP=-;
MATRIX
'''.format(aln_info[0], aln_info[1]))
        
        for rec in contents:
            fh.write("\n{} {}".format(rec.description, str(rec.seq)))
        fh.write("\n;\nEnd;")
        
    shutil.move(out, ndir)
                       
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
    ndir = os.path.join(curpath, "Nexus-Files")
    if not os.path.exists(ndir):
        os.mkdir(ndir)
        
    pdir = os.path.join(curpath, "Phylip-Files")
    if not os.path.exists(pdir):
        os.mkdir(pdir)

    return ndir, pdir

def main():
    tb = datetime.now()
    args = get_args()
    
    ndir, pdir = make_dirs(args.outdir)
    
    os.chdir(args.indir)
    flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
    
    print("\n\nFound {} fasta files to convert.\n\n".format(len(flist)))
    
    for f in flist:
        prefix = f.split('.')[0]
        contents, aln_info = read_fasta(f)
        convert_to_phy(contents, aln_info, prefix, pdir)
        convert_to_nex(contents, aln_info, prefix, ndir)

    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
if __name__ == '__main__':
    main()


