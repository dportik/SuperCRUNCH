''''
SuperCRUNCH: Fasta_Convert module

Usage: python Fasta_Convert.py -i [directory with all fasta files]

    Fasta_Convert: Converts a directory of aligned fasta files into both phylip and nexus formats.
    Moves new files into corresponding output directories. Note this should only be used after
    records have been renamed with species or accession names, as the original NCBI description
    lines will cause problems (too long and with many spaces or '|' within).
    
    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.


    Output_Phylip_Files:
           [NAME].phy - The phylip format alignment.
                                           
    Output_Nexus_Files:
           [NAME].nex - The nexus format alignment
           
-------------------------
For Python 2.7
Python modules required:
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
from Bio import SeqIO

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Fasta_Convert: Converts a directory of aligned fasta files into both phylip and nexus formats.
    Moves new files into corresponding output directories. Note this should only be used after
    records have been renamed with species or accession names, as the original NCBI description
    lines will cause problems (too long and with many spaces or '|' within).
    
    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.


    Output_Phylip_Files:
           [NAME].phy - The phylip format alignment.
                                           
    Output_Nexus_Files:
           [NAME].nex - The nexus format alignment
        DEPENDENCIES: Python: BioPython.
           	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the fasta format alignments to convert. Follow labeling format: NAME.fasta")
    return parser.parse_args()

def make_nested_dir(in_dir, string):
    dirname = in_dir+'/'+string
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    return dirname

def split_name(string, index, delimiter):
    index = int(index)
    name = string.split(delimiter)[index]
    return name

def read_fasta(f):
    contents = list(SeqIO.parse(f,'fasta'))
    aln_info = [len(contents), len(contents[0].seq)]            
    return contents, aln_info

def convert_to_phy(contents, aln_info, prefix, out_dir):
    print "\tConverting {}.fasta into phylip format".format(prefix)
    out = "{}.phy".format(prefix)
    with open(out, 'a') as fh_out:
        fh_out.write("{0} {1}".format(aln_info[0],aln_info[1]))
        for rec in contents:
            fh_out.write("\n{} {}".format(rec.description, str(rec.seq)))
    shutil.move(out, out_dir)
    
def convert_to_nex(contents, aln_info, prefix, out_dir):
    print "\tConverting {}.fasta into nexus format\n".format(prefix)
    out = "{}.nex".format(prefix)
    with open(out, 'a') as fh_out:
        fh_out.write('''#NEXUS 
BEGIN DATA;
	DIMENSIONS  NTAX={0} NCHAR={1};
	FORMAT DATATYPE=DNA  MISSING=N GAP=-;
MATRIX
'''.format(aln_info[0],aln_info[1]))
        for rec in contents:
            fh_out.write("\n{} {}".format(rec.description, str(rec.seq)))
        fh_out.write("\n;\nEnd;")
    shutil.move(out, out_dir)
                       
#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    os.chdir(args.in_dir)
    fasta_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])
    phy_dir = make_nested_dir(args.in_dir, 'Output_Phylip_Files')
    nex_dir = make_nested_dir(args.in_dir, 'Output_Nexus_Files')
    print "\n\nFound {} fasta files to convert.\n\n".format(len(fasta_list))
    for f in fasta_list:
        contents, aln_info = read_fasta(f)
        prefix = split_name(f, 0, '.')
        convert_to_phy(contents, aln_info, prefix, phy_dir)
        convert_to_nex(contents, aln_info, prefix, nex_dir)
    
if __name__ == '__main__':
    main()


