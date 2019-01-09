''''
SuperCRUNCH: Trim_Alignments module

Usage: python Trim_Alignments.py  -i [directory with all alignment files] (REQUIRED)
                                  -a [trimal method (choices: gt, noallgaps, both)] (REQUIRED)
                                  -f [output format (choices: fasta, nexus, or phylip)] (REQUIRED)
                                  --gt [value for gap threshold arg in trimal (0.01-1)] (Optional)

	Trim_Alignments: Use the program trimal to batch trim all alignment files in a directory.
    There are three options (-a) for using trimal in this script. The first is to use the 'gt' method. The
    default value is set to 0.05, meaning for a column to be kept 95% of the
    sequences must not have a gap. This value can be changed using the --gt
    flag with your own value selected. The second option is 'noallgaps', which removes
    any columns composed entirely of gaps. The third option 'both' runs the gt method, followed
    by the noallgaps method.
    
    Input alignment file formats are auto-detected and can be include fasta, nexus, or 
    phylip formats. The input alignment files must be labeled with one of the following 
    extensions to be read: NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, NAME.phylip, or NAME.phy.
    
    The output format must be specified using the -f argument (choices: fasta, nexus, or phylip). 

    Timmed alignment output files will be moved to following output directory
    created in the main alignment directory:

    /Output_Trimmed_Alignments:
           [NAME]_trimmed.[phy, fasta, or nex extension]
    
-------------------------
For Python 2.7
Dependencies:
	-trimal (installed in path)
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
import subprocess as sp
import shutil
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
	Trim_Alignments: Use the program trimal to batch trim all alignment files in a directory.
    There are three options (-a) for using trimal in this script. The first is to use the 'gt' method. The
    default value is set to 0.05, meaning for a column to be kept 95% of the
    sequences must not have a gap. This value can be changed using the --gt
    flag with your own value selected. The second option is 'noallgaps', which removes
    any columns composed entirely of gaps. The third option 'both' runs the gt method, followed
    by the noallgaps method.
    
    Input alignment file formats are auto-detected and can be include fasta, nexus, or 
    phylip formats. The input alignment files must be labeled with one of the following 
    extensions to be read: NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, NAME.phylip, or NAME.phy.
    
    The output format must be specified using the -f argument (choices: fasta, nexus, or phylip). 

    Timmed alignment output files will be moved to following output directory
    created in the main alignment directory:

    /Output_Trimmed_Alignments:
           [NAME]_trimmed.[phy, fasta, or nex extension]

    DEPENDENCIES: Executables in path: trimal (with name identical to that listed here).
        ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files. File formats are auto-detected and can be include fasta, nexus, or phylip formats, but input alignment files must be labeled with one of the following extensions: NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, NAME.phylip, or NAME.phy.")
    parser.add_argument("-f", "--format", required=True, choices=["fasta","nexus","phylip"], help="REQUIRED: Specify the output file format for trimmed alignment.")
    parser.add_argument("-a", "--analysis", required=True, choices=["gt","noallgaps","both"], help="REQUIRED: Specify the trimal method for trimming alignments.")
    parser.add_argument("--gt", default="0.05", help="OPTIONAL: Specify the gt (gap threshold) value for trimal - the minimum fraction of sequences without a gap. (default=0.05)")
    return parser.parse_args()

def trim_align_gt(aln_file, out_format, gt):
    print "\n\nTrimming alignment {}:".format(aln_file)
    prefix = aln_file.split('.')[0]
    if out_format == "fasta":
        call_string = "trimal -in {0} -out {1}_trimmed.fasta -fasta -gt {2}".format(aln_file, prefix, gt)
    elif out_format == "nexus":
        call_string = "trimal -in {0} -out {1}_trimmed.nex -nexus -gt {2}".format(aln_file, prefix, gt)
    elif out_format == "phylip":
        call_string = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -gt {2}".format(aln_file, prefix, gt)
    print call_string
    proc = sp.call(call_string, shell=True)
    
def trim_align_nag(aln_file, out_format):
    print "\n\nTrimming alignment {}:".format(aln_file)
    prefix = aln_file.split('.')[0]
    if out_format == "fasta":
        call_string = "trimal -in {0} -out {1}_trimmed.fasta -fasta -noallgaps".format(aln_file, prefix)
    elif out_format == "nexus":
        call_string = "trimal -in {0} -out {1}_trimmed.nex -nexus -noallgaps".format(aln_file, prefix)
    elif out_format == "phylip":
        call_string = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -noallgaps".format(aln_file, prefix)
    print call_string
    proc = sp.call(call_string, shell=True)


def trim_align_both(aln_file, out_format, gt):
    print "\n\nTrimming alignment {}:".format(aln_file)
    prefix = aln_file.split('.')[0]
    if out_format == "fasta":
        call_string1 = "trimal -in {0} -out {1}_temp.fasta -fasta -gt {2}".format(aln_file, prefix, gt)
    elif out_format == "nexus":
        call_string1 = "trimal -in {0} -out {1}_temp.nex -nexus -gt {2}".format(aln_file, prefix, gt)
    elif out_format == "phylip":
        call_string1 = "trimal -in {0} -out {1}_temp.phy -phylip_paml -gt {2}".format(aln_file, prefix, gt)
    print call_string1
    proc = sp.call(call_string1, shell=True)

    if out_format == "fasta":
        call_string2 = "trimal -in {0} -out {1}_trimmed.fasta -fasta -noallgaps".format("{0}_temp.fasta".format(prefix), prefix)
    elif out_format == "nexus":
        call_string2 = "trimal -in {0} -out {1}_trimmed.nex -nexus -noallgaps".format("{0}_temp.nexus".format(prefix), prefix)
    elif out_format == "phylip":
        call_string2 = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -noallgaps".format("{0}_temp.phy".format(prefix), prefix)
    print call_string2
    proc = sp.call(call_string2, shell=True)

    temp_list = [t for t in os.listdir('.') if t.endswith('_temp.fasta') or t.endswith('_temp.nex') or t.endswith('_temp.phy')]
    for t in temp_list:
        os.remove(t)

def main():
    tb = datetime.now()
    args = get_args()
    os.chdir(args.in_dir)
    out_dir = "Output_Trimmed_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    f_list = sorted([f for f in os.listdir('.') if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.nexus') or f.endswith('.nex') or f.endswith('.phy') or f.endswith('.phylip')])
    print "\nBeginning trimming of alignments.\n"
    for f in f_list:
        if args.analysis == "gt":
            trim_align_gt(f, args.format, args.gt)
        elif args.analysis == "noallgaps":
            trim_align_nag(f, args.format)
        elif args.analysis == "both":
            trim_align_both(f, args.format, args.gt)
        out_list = [o for o in os.listdir('.') if o.endswith('_trimmed.fasta') or o.endswith('_trimmed.nex') or o.endswith('_trimmed.phy')]
        for out in out_list:
            shutil.move(out, out_dir)
    tf = datetime.now()
    te = tf - tb
    print "\n\n\nFinished trimming alignments. Elapsed time: {} (H:M:S)\n\n".format(te)


if __name__ == '__main__':
    main()
