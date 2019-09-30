''''
SuperCRUNCH: Trim_Alignments_Trimal module

	Trim_Alignments_Trimal: Use the program trimal to batch trim all alignment files in a 
    directory. There are four options (-a) for using trimal in this script. The 
    first is to use the 'gt' method. The default value is set to 0.05, meaning for 
    a column to be kept 95% of the sequences must not contain a gap. This value can 
    be changed using the --gt flag with your own value selected. The second option 
    is 'noallgaps', which removes any columns composed entirely of gaps. The third 
    option 'both' runs the gt method, followed by the noallgaps method. The fourth 
    option is to run the 'gappyout' method.
    
    Input alignment file formats are auto-detected and can be include fasta, nexus, 
    or phylip formats. The input alignment files must be labeled with one of the 
    following extensions to be read: NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, 
    NAME.phylip, or NAME.phy. *NOTE: these files should have the description lines
    relabeled prior to trimming, otherwise only the accession number will be written
    to the new trimmed files.
    
    The output format must be specified using the -f argument (choices: fasta, 
    nexus, or phylip). 

    Trimmed alignment output files will be moved to following output directory specified.

-------------------------
Compatible with Python 2.7 & 3.7
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

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
	Trim_Alignments_Trimal: Use the program trimal to batch trim all alignment files in a 
    directory. There are four options (-a) for using trimal in this script. The 
    first is to use the 'gt' method. The default value is set to 0.05, meaning for 
    a column to be kept 95% of the sequences must not contain a gap. This value can 
    be changed using the --gt flag with your own value selected. The second option 
    is 'noallgaps', which removes any columns composed entirely of gaps. The third 
    option 'both' runs the gt method, followed by the noallgaps method. The fourth 
    option is to run the 'gappyout' method.
    
    Input alignment file formats are auto-detected and can be include fasta, nexus, 
    or phylip formats. The input alignment files must be labeled with one of the 
    following extensions to be read: NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, 
    NAME.phylip, or NAME.phy. *NOTE: these files should have the description lines
    relabeled prior to trimming, otherwise only the accession number will be written
    to the new trimmed files.
    
    The output format must be specified using the -f argument (choices: fasta, 
    nexus, or phylip). 

    Trimmed alignment output files will be moved to following output directory specified.

    DEPENDENCIES: Executables in path: trimal.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input alignment files. File formats are auto-detected and "
                            "can be include fasta, nexus, or phylip formats, but input alignment "
                            "files must be labeled with one of the following extensions: "
                            "NAME.fasta, NAME.fa, NAME.nexus, NAME.nex, NAME.phylip, or NAME.phy.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-f", "--format",
                            required=True,
                            choices=["fasta","nexus","phylip"],
                            help="REQUIRED: Specify the output file format for trimmed alignment.")
    
    parser.add_argument("-a", "--analysis",
                            required=True,
                            choices=["gt","noallgaps","both","gappyout"],
                            help="REQUIRED: Specify the trimal method for trimming alignments.")
    
    parser.add_argument("--gt",
                            default="0.05",
                            help="OPTIONAL: Specify the gt (gap threshold) value for trimal - "
                            "the minimum fraction of sequences without a gap. (default=0.05)")
    
    return parser.parse_args()

def trim_align_gt(aln_file, out_format, gt):
    """
    Trim alns in correct format using 'gap threshold' option 
    and unless user provides value, use default of 0.05.
    """
    print("\n\nTrimming alignment {}:".format(aln_file))
    prefix = aln_file.split('.')[0]
    
    if out_format == "fasta":
        cmd = "trimal -in {0} -out {1}_trimmed.fasta -fasta -gt {2}".format(aln_file, prefix, gt)
        
    elif out_format == "nexus":
        cmd = "trimal -in {0} -out {1}_trimmed.nex -nexus -gt {2}".format(aln_file, prefix, gt)
        
    elif out_format == "phylip":
        cmd = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -gt {2}".format(aln_file, prefix, gt)
        
    print(cmd)
    proc = sp.call(cmd, shell=True)
    
def trim_align_nag(aln_file, out_format):
    """
    Trim alns in correct format using 'noallgaps' option.
    """
    print("\n\nTrimming alignment {}:".format(aln_file))
    prefix = aln_file.split('.')[0]
    
    if out_format == "fasta":
        cmd = "trimal -in {0} -out {1}_trimmed.fasta -fasta -noallgaps".format(aln_file, prefix)
        
    elif out_format == "nexus":
        cmd = "trimal -in {0} -out {1}_trimmed.nex -nexus -noallgaps".format(aln_file, prefix)
        
    elif out_format == "phylip":
        cmd = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -noallgaps".format(aln_file, prefix)
        
    print(cmd)
    proc = sp.call(cmd, shell=True)

def trim_align_go(aln_file, out_format):
    """
    Trim alns in correct format using 'gappyout' option.
    """
    print("\n\nTrimming alignment {}:".format(aln_file))
    prefix = aln_file.split('.')[0]
    
    if out_format == "fasta":
        cmd = "trimal -in {0} -out {1}_trimmed.fasta -fasta -gappyout".format(aln_file, prefix)
        
    elif out_format == "nexus":
        cmd = "trimal -in {0} -out {1}_trimmed.nex -nexus -gappyout".format(aln_file, prefix)
        
    elif out_format == "phylip":
        cmd = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -gappyout".format(aln_file, prefix)
        
    print(cmd)
    proc = sp.call(cmd, shell=True)

def trim_align_both(aln_file, out_format, gt):
    """
    Trim alns in correct format using gap threshold option 
    and unless user provides value, use default of 0.05. 
    Follow with noallgaps option to be sure no gapped columns.
    """
    print("\n\nTrimming alignment {}:".format(aln_file))
    prefix = aln_file.split('.')[0]
    
    if out_format == "fasta":
        cmd1 = "trimal -in {0} -out {1}_temp.fasta -fasta -gt {2}".format(aln_file, prefix, gt)
        
    elif out_format == "nexus":
        cmd1 = "trimal -in {0} -out {1}_temp.nex -nexus -gt {2}".format(aln_file, prefix, gt)
        
    elif out_format == "phylip":
        cmd1 = "trimal -in {0} -out {1}_temp.phy -phylip_paml -gt {2}".format(aln_file, prefix, gt)
        
    print(cmd1)
    proc = sp.call(cmd1, shell=True)

    if out_format == "fasta":
        cmd2 = "trimal -in {0} -out {1}_trimmed.fasta -fasta -noallgaps".format("{0}_temp.fasta".format(prefix), prefix)
        
    elif out_format == "nexus":
        cmd2 = "trimal -in {0} -out {1}_trimmed.nex -nexus -noallgaps".format("{0}_temp.nexus".format(prefix), prefix)
        
    elif out_format == "phylip":
        cmd2 = "trimal -in {0} -out {1}_trimmed.phy -phylip_paml -noallgaps".format("{0}_temp.phy".format(prefix), prefix)
        
    print(cmd2)
    proc = sp.call(cmd2, shell=True)

    [os.remove(t) for t in os.listdir('.') if t.endswith(('_temp.fasta', '_temp.nex', '_temp.phy'))]
    
def make_dirs(outdir, f_format):
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    os.chdir(outdir)
    curpath = os.getcwd()
    
    tdir = os.path.join(curpath, "Trimmed-{}-Files".format(f_format))
    if not os.path.exists(tdir):
        os.mkdir(tdir)

    return tdir
    
def main():
    tb = datetime.now()
    args = get_args()

    tdir = make_dirs(args.outdir, args.format)
    
    os.chdir(args.indir)
            
    flist = sorted([f for f in os.listdir('.')
                         if f.endswith(('.fasta', '.fa',
                                            '.nexus', '.nex',
                                            '.phy', '.phylip'))])
    
    print("\nBeginning trimming of alignments.\n")
    
    for f in flist:
        if args.analysis == "gt":
            trim_align_gt(f, args.format, args.gt)
            
        elif args.analysis == "noallgaps":
            trim_align_nag(f, args.format)
            
        elif args.analysis == "both":
            trim_align_both(f, args.format, args.gt)
            
        elif args.analysis == "gappyout":
            trim_align_go(f, args.format)
            
        [shutil.move(o, tdir) for o in os.listdir('.')
             if o.endswith(('trimmed.fasta', 'trimmed.nex', 'trimmed.phy'))]
        
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished trimming alignments. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    


if __name__ == '__main__':
    main()
