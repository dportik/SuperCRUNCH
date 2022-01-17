''''
SuperCRUNCH: Align module

    Align: Perform alignments for a directory of unaligned fasta files using 
    MAFFT, MACSE, MUSCLE, or Clustal-O. 

    Alignment method is specified using the -a flag (options: mafft, macse, 
    muscle, clustalo). To run mafft, muscle, and clustalo the executables must
    be in path with names matching those here. To run macse, the user must provide
    the full path to a macse.jar file using the --mpath flag. Selecting 'all' 
    will run mafft, muscle, and clustalo (but not macse) sequentially.

    Several of the aligment methods will either run under the default 
    settings (muscle) or the auto select settings (mafft, clustalo). 
    The optional flag --accurate can be used for mafft (invokes FFT-NS-i; 
    e.g., --retree 2 --maxiterate 1000) and clustalo (enables --iter=5, 
    in which the guide-tree and HMM each undergo 5 iterations, rather
    than only one).
    
    The translation alignment method of macse requires a translation table, 
    which here defaults to standard code unless specified using the --table 
    flag. Many (but not all) NCBI translation table options are available, and 
    can be selected using integers or the shortcut terms provided. The --pass_fail
    can be used to signal that potentially two fasta files should be dual aligned
    for a given locus (see documentation for details). Additional 
    memory (in GB) can be assigned to macse using the --mem flag. The 
    --accurate flag can be used for macse 2.0+, and adds the commands 
    -local_realign_init 0.9 and -local_realign_dec 0.9, which will increase 
    alignment time but improve accuracy, sometimes considerably. 
 
    Output files vary between aligners but will be moved to an aligner-specific
    output directory in the main output directory specified by the user.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
Dependencies:
	-mafft (in path)
    -macse (requires jar file)
    -muscle (in path)
    -clustalo (in path)
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
 2019
Distributed under the 
GNU General Public Lincense
'''

import argparse
import os
import subprocess as sp
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Align: Perform alignments for a directory of unaligned fasta files using 
    MAFFT, MACSE, MUSCLE, or Clustal-O. 

    Alignment method is specified using the -a flag (options: mafft, macse, 
    muscle, clustalo). To run mafft, muscle, and clustalo the executables must
    be in path with names matching those here. To run macse, the user must provide
    the full path to a macse.jar file using the --mpath flag. Selecting 'all' 
    will run mafft, muscle, and clustalo (but not macse) sequentially.

    Several of the aligment methods will either run under the default 
    settings (muscle) or the auto select settings (mafft, clustalo). 
    The optional flag --accurate can be used for mafft (invokes FFT-NS-i; 
    e.g., --retree 2 --maxiterate 1000) and clustalo (enables --iter=5, 
    in which the guide-tree and HMM each undergo 5 iterations, rather
    than only one).
    
    The translation alignment method of macse requires a translation table, 
    which here defaults to standard code unless specified using the --table 
    flag. Many (but not all) NCBI translation table options are available, and 
    can be selected using integers or the shortcut terms provided. The --pass_fail
    can be used to signal that potentially two fasta files should be dual aligned
    for a given locus (see documentation for details). Additional 
    memory (in GB) can be assigned to macse using the --mem flag. The 
    --accurate flag can be used for macse 2.0+, and adds the commands 
    -local_realign_init 0.9 and -local_realign_dec 0.9, which will increase 
    alignment time but improve accuracy, sometimes considerably. 
 
    Output files vary between aligners but will be moved to an aligner-specific
    output directory in the main output directory specified by the user.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
    
    DEPENDENCIES: Python: BioPython; Executables in path: mafft, muscle, 
    clustalo; macse jar file required.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-a", "--aln",
                            required=True,
                            choices=["mafft","macse","muscle", "clustalo","all"],
                            help="REQUIRED: Specify whether alignment is by mafft, macse, "
                            "muscle, or clustalo. If macse must provide flags --mpath and "
                            "--table. Selecting 'all' will run mafft, muscle, and clustalo "
                            "sequentially.")
    
    parser.add_argument("--accurate",
                            action='store_true',
                            help="OPTIONAL: Specifies mafft, clustalo, or macse to use "
                            "more thorough search settings.")
    
    parser.add_argument("--threads",
                            default=None,
                            help="OPTIONAL: Specifies number of threads to use in mafft "
                            "and/or clustalo.")
    
    parser.add_argument("--mpath",
                            default=None,
                            help="MACSE REQUIRED: Full path to a macse.jar file.")
    
    parser.add_argument("--table",
                            default="standard",
                            choices=["standard", "vertmtdna", "invertmtdna", "yeastmtdna",
                                         "plastid", "1", "2", "3", "4", "5", "6", "9",
                                         "10", "11", "12", "13", "14", "15", "16", "21",
                                         "22", "23"],
                            help="MACSE REQUIRED: Specifies translation table for macse.")
    
    parser.add_argument("--mem",
                            default=None,
                            help="MACSE OPTIONAL: An integer to assign additional memory "
                            "to macse (in GB), default=1.")
    
    parser.add_argument("--pass_fail",
                            action='store_true',
                            help="MACSE OPTIONAL: Specifies macse to use two fasta files "
                            "for dual file alignment. See documentation for details.")
    
    
    return parser.parse_args()

def table_dict(symbol):
    """
    Get MACSE translation table value from user choice.
    Uses dictionaries to look up string that will be
    used in the macse call.
    """
    table_dictv1 = {"standard":"1", "vertmtdna":"2", "invertmtdna":"5", "yeastmtdna":"3",
                        "plastid":"11", "1":"1", "2":"2", "3":"3", "4":"4", "5":"5", "6":"6",
                        "9":"9", "10":"10", "11":"11", "12":"12", "13":"13", "14":"14",
                        "15":"15", "16":"16", "21":"21", "22":"22", "23":"23"}
        
    table_dictv2 = {"standard":"01_The_Standard_Code",
                        "vertmtdna":"02_The_Vertebrate_Mitochondrial_Code",
                        "invertmtdna":"05_The_Invertebrate_Mitochondrial_Code",
                        "yeastmtdna":"03_The_Yeast_Mitochondrial_Code",
                        "plastid":"11_The_Bacterial_Archaeal_and_Plant_Plastid_Code",
                        "1":"01_The_Standard_Code",
                        "2":"02_The_Vertebrate_Mitochondrial_Code",
                        "3":"03_The_Yeast_Mitochondrial_Code",
                        "4":"04_The_Mold_Protozoan_and_Coelenterate_Mitochondrial_Code_and_the_Mycoplasma_Spiroplasma_Code",
                        "5":"05_The_Invertebrate_Mitochondrial_Code",
                        "6":"06_The_Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code",
                        "9":"09_The_Echinoderm_and_Flatworm_Mitochondrial_Code",
                        "10":"10_The_Euplotid_Nuclear_Code",
                        "11":"11_The_Bacterial_Archaeal_and_Plant_Plastid_Code",
                        "12":"12_The_Alternative_Yeast_Nuclear_Code",
                        "13":"13_The_Ascidian_Mitochondrial_Code",
                        "14":"14_The_Alternative_Flatworm_Mitochondrial_Code",
                        "15":"15_Blepharisma_Nuclear_Code",
                        "16":"16_Chlorophycean_Mitochondrial_Code",
                        "21":"21_Trematode_Mitochondrial_Code",
                        "22":"22_Scenedesmus_obliquus_mitochondrial_Code",
                        "23":"23_Thraustochytrium_Mitochondrial_Code"}
        
    table_val = table_dictv1[symbol]
    
    return table_val


def pass_fail_finder(flist):
    """
    Finds files to use for the paired or unpaired
    macse alignments. Paired analyses will have two files
    which end in '_Passed.fasta' and '_Failed.fasta'. The
    unpaired analyses will only have a '_Passed.fasta' file,
    and any loci with only failed sequences ('_Failed.fasta')
    are excluded from alignments. Uses set methods to identify
    the paired and unpaired file names after searching for
    the files, then returns two lists of the prefixes for
    the paired and unpaired analyses.
    """
    #first get prefixes from sets of files
    prefixes1 = set([f.split('_Passed')[0] for f in flist if f.endswith('_Passed.fasta')])
    prefixes2 = set([f.split('_Failed')[0] for f in flist if f.endswith('_Failed.fasta')])

    #create sorted list of prefixes that only
    #have a '_Passed.fasta' file present
    pass_only = sorted(prefixes1 - prefixes2)
    
    #create sorted list of prefixes that only
    #have a both a '_Passed.fasta' and a
    #'_Failed.fasta' file present
    pass_fail = sorted(prefixes1 & prefixes2)

    #create sorted list of prefixes that only
    #have a '_Failed.fasta' file present
    fail_only = sorted(prefixes2 - prefixes1)
    
    if pass_fail:
        print("\nFound {} loci with pass and fail sequence pairs:".format(len(pass_fail)))
        for p in pass_fail:
            print("\t{}".format(p))
    
    if pass_only:
        print("\nFound {} loci with only pass sequences:".format(len(pass_only)))
        for p in pass_only:
            print("\t{}".format(p))
            
    if fail_only:
        print("\nFound {} loci with only fail sequences, which will be excluded from alignment:"
                  .format(len(fail_only)))
        for f in fail_only:
            print("\t{}".format(f))
    
    return pass_only, pass_fail


def get_cmd(f, aln, accurate, threads, mpath, table, mem):
    """
    Generate relevant call string for file f based on
    arguments supplied. Returns a call string for f. This
    is used to generate a command for all methods except
    for the macse --pass_fail method. 
    """
    #params used by multiple alignment methods
    #set threads to 1 if not set
    if threads is None:
        threads = 1
    #get a prefix to label output files
    #based on input file name (split by period)
    prefix = f.split('.')[0]
    
    if aln == "mafft":
        outname = "{}_mafft_temp.fasta".format(prefix)
        if accurate is False:
            cmd = ("mafft --thread {0} --auto {1} > {2}"
                               .format(threads, f, outname))
        elif accurate is True:
            cmd = ("mafft --thread {0} --retree 2 --maxiterate 1000 {1} > {2}"
                               .format(threads, f, outname))
                
    elif aln == "muscle":
        outname = "{}_muscle_temp.fasta".format(prefix)
        cmd = "muscle -in {0} -out {1}".format(f, outname)
        

    elif aln == "clustalo":
        outname = "{}_clustalo_temp.fasta".format(prefix)
        if accurate is False:
            cmd = ("clustalo -i {0} -o {1} --auto -v --threads={2} --output-order=tree-order --force"
                       .format(f, outname, threads))
        elif accurate is True:
            cmd = ("clustalo -i {0} -o {1} --full --full-iter --iter=5 -v --threads={2} --cluster-size=500 --output-order=tree-order --force"
                       .format(f, outname, threads))

    elif aln == "macse":
        #macse specific parameters
        #assign memory of 1GB if none set
        if mem is None:
            mem = 1
        #get translation table string (default is standard)
        tcode = table_dict(table)
        if accurate is False:
            cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} "
                               .format(mem, mpath, tcode, f))
        elif accurate is True:
            cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -local_realign_init 0.9 -local_realign_dec 0.9 "
                               .format(mem, mpath, tcode, f))

    return cmd
                
def get_cmds_macse_pass_fail(flist, accurate, mpath, table, mem):
    """
    Generate a list of commands for macse --pass_fail
    alignments, based on the file list and arguments
    supplied. Relies on the pass_fail_finder() function
    to identify paired and/or unpaired macse files and
    the produces the proper call strings.
    """
    #macse specific parameters
    #assign memory of 1GB if none set
    if mem is None:
        mem = 1
    #get translation table string (default is standard)
    tcode = table_dict(table)

    #find out which files have pass fail pairs,
    #and pass only files (ignores fail only)
    #returns sorted lists of prefixes
    pass_only, pass_fail = pass_fail_finder(flist)

    #set up list to populate with commands
    commands = []

    #files with pass only may not exist, check first
    if pass_only:
        for prefix in pass_only:
            #re-constitute original filename
            f = '{}_Passed.fasta'.format(prefix)

            if accurate is False:
                cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} "
                                   .format(mem, mpath, tcode, f))
                commands.append(cmd)
                
            elif accurate is True:
                cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -local_realign_init 0.9 -local_realign_dec 0.9 "
                                   .format(mem, mpath, tcode, f))
                commands.append(cmd)
                
    #files with pass fail pairs may not exist, check first
    if pass_fail:
        for prefix in pass_fail:
            #re-constitute original filenames
            fpass = '{}_Passed.fasta'.format(prefix)
            ffail = '{}_Failed.fasta'.format(prefix)

            if accurate is False:
                cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -seq_lr {4} "
                           .format(mem, mpath, tcode, fpass, ffail))
                commands.append(cmd)
                
            elif accurate is True:
                cmd = ("java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -seq_lr {4} -local_realign_init 0.9 -local_realign_dec 0.9 "
                           .format(mem, mpath, tcode, fpass, ffail))
                commands.append(cmd)

    return commands
        
def get_all_commands(flist, aln, accurate, threads, mpath, table, mem, pass_fail):
    """
    Generates a complete set of commands to run alignments
    for all files in flist. In all cases, returns a list of 
    strings which will run the alignment method for a given 
    file using the sp.call() method with the shell argument.
    """
    
    if aln == "macse" and pass_fail is True:
        commands = get_cmds_macse_pass_fail(flist, accurate, mpath, table, mem)

    elif aln == "all":
        commands = []
        commands.extend([get_cmd(f, "mafft", accurate, threads, mpath, table, mem) for f in flist])
        commands.extend([get_cmd(f, "clustalo", accurate, threads, mpath, table, mem) for f in flist])
        commands.extend([get_cmd(f, "muscle", accurate, threads, mpath, table, mem) for f in flist])

    else:
        commands = [get_cmd(f, aln, accurate, threads, mpath, table, mem) for f in flist]

    return commands
        

def reformat(f, aln, outdir):
    """
    Many of the output fasta files from the alignment methods
    are on the messy side. Here, the temporary output alignment
    file is read using biopython and re-written in a consistent
    fasta format. The output file is named using the prefix of 
    the original fasta file and the alignment method chosen. 
    The new alignment file is moved to the proper output directory
    and the temporary file is removed.
    """
    fdict = SeqIO.index(f, "fasta")
    
    prefix = f.split("_{}_".format(aln))[0]
    outname = "{0}_{1}_Aligned.fasta".format(prefix, aln.upper())
    
    with open(outname, 'a') as fh:
        for record in fdict:
            newseq = fdict[record].seq.upper()
            fh.write( ">{}\n{}\n".format(fdict[record].description, newseq))
            
    os.remove(f)
    
    try:
        shutil.move(outname, outdir)
    except:
        os.remove(os.path.join(outdir, outname))
        shutil.move(outname, outdir)
        print("\n\nWARNING: File '{}' already exists in:\n\t{}.\n\tReplacing previous version.".format(outname, outdir))

def macse_reformat(f, outdir):
    """
    In cases of orf errors, macse inserts a ! character
    which is 'illegal' in most downstream analyses. This
    function opens the temporary output alignment file and
    replaces any ! characters with an N instead, then writes
    the contents to a new alignment file named using the 
    prefix of the original fasta file and the alignment method.
    The new alignment file is moved to the proper output directory
    and the original alignment files (amino-acid translated and
    nucleotid) are moved to their own output directory.
    """
    with open(f, 'r') as fh:
        cleaned_lines = [line.replace('!','N') for line in fh]
        
    prefix = f.split('_Passed_NT')[0]
    outname = "{}_MACSE_Aligned.fasta".format(prefix)
    
    with open(outname, 'a') as fh:
        for l in cleaned_lines:
            fh.write(l)
            
    try:
        shutil.move(outname, outdir)
    except:
        os.remove(os.path.join(outdir, outname))
        shutil.move(outname, outdir)
        print("\n\nWARNING: File '{}' already exists in:\n\t{}.\n\tReplacing previous version.".format(outname, outdir))

def cleanup(aln, mafftdir, muscledir, clustaldir, macsecdir, macseodir, logpath):
    """
    Based on alignment method selected by user (aln), find the temporary
    output alignment file and produce the final output file using either
    the reformat() or macse_reformat() functions. 
    """
    if aln == "mafft":
        temps = [f for f in os.listdir('.') if f.endswith("_mafft_temp.fasta") and os.stat(f).st_size != 0]
        if not temps:
            message = ("\n\n\nERROR:\tNo mafft output files were produced. Please ensure mafft is installed in path.\n\n\n")
            write_log(logpath, message)
            raise ValueError(message)
        else:
            for f in temps:
                reformat(f, aln, mafftdir) 
        
    elif aln == "muscle":
        temps = [f for f in os.listdir('.') if f.endswith("_muscle_temp.fasta")]
        if not temps:
            message = ("\n\n\nERROR:\tNo muscle output files were produced. Please ensure muscle is installed in path.\n\n\n")
            write_log(logpath, message)
            raise ValueError(message)
        else:
            for f in temps:
                reformat(f, aln, muscledir)
                
    elif aln == "clustalo":
        temps = [f for f in os.listdir('.') if f.endswith("_clustalo_temp.fasta")]
        if not temps:
            message = ("\n\n\nERROR:\tNo clustalo output files were produced. Please ensure clustalo is installed in path.\n\n\n")
            write_log(logpath, message)
            raise ValueError(message)
        else:
            for f in temps:
                reformat(f, aln, clustaldir)
        
    elif aln == "macse":
        temps = [f for f in os.listdir('.') if f.endswith("_NT.fasta")]
        if not temps:
            message = ("\n\n\nERROR:\tNo macse output files were produced. Please ensure macse jar path is correct.\n\n\n")
            write_log(logpath, message)
            raise ValueError(message)
        else:
            for f in temps:
                macse_reformat(f, macsecdir)
        for f in [x for x in os.listdir('.') if x.endswith(("_AA.fasta", "_NT.fasta"))]:
            try:
                shutil.move(f, macseodir)
            except:
                os.remove(os.path.join(macseodir, f))
                shutil.move(f, macseodir)
                print("\n\nWARNING: File '{}' already exists in:\n\t{}.\n\tReplacing previous version.".format(f, macseodir))
        
    elif aln == "all":
        temps1 = [f for f in os.listdir('.') if f.endswith("_mafft_temp.fasta") and os.stat(f).st_size != 0]
        for f in temps1:
            reformat(f, "mafft", mafftdir)
                
        temps2 = [f for f in os.listdir('.') if f.endswith("_muscle_temp.fasta")]
        for f in temps2:
            reformat(f, "muscle", muscledir)
                
        temps3 = [f for f in os.listdir('.') if f.endswith("_clustalo_temp.fasta")]
        for f in temps3:
            reformat(f, "clustalo", clustaldir)

def clear_dir(p, logpath):
    message = ("WARNING: An output directory already exists:\n\t'{}'\n\tRemoving "
                   "all files in this directory to prevent errors.\n".format(p))
    print("\n{}".format(message))
    write_log(logpath, message)
    shutil.rmtree(p)
    os.mkdir(p)

def make_dirs(outdir, aln, logpath):
    """
    Creates path names for directories for all alignment
    methods, but only creates the relevant directories
    based on the aln argument supplied by the user. Returns
    all directory paths (even those not created) to allow
    output files to be moved to the proper directories 
    during cleanup steps.
    """
    #get current path
    os.chdir(outdir)
    curpath = os.getcwd()
    
    #create paths using os.path.join() to avoid any issues
    mafftdir = os.path.join(curpath, "Alignments-MAFFT")
    muscledir = os.path.join(curpath, "Alignments-MUSCLE")
    clustaldir = os.path.join(curpath, "Alignments-CLUSTALO")
    macsedir = os.path.join(curpath, "Alignments-MACSE")    
    macsecdir = os.path.join(curpath, "Alignments-MACSE", "Cleaned-Alignments")
    macseodir = os.path.join(curpath, "Alignments-MACSE", "Additional-Outputs")

    #create the directories using the above paths
    #if the alignment method was selected by user
    if aln == "mafft":
        if not os.path.exists(mafftdir):
            os.mkdir(mafftdir)
            
    elif aln == "muscle":
        if not os.path.exists(muscledir):
            os.mkdir(muscledir)
            
    elif aln == "clustalo":
        if not os.path.exists(clustaldir):
            os.mkdir(clustaldir)

    elif aln == "macse":
        for p in [macsedir, macsecdir, macseodir]:
            if not os.path.exists(p):
                os.mkdir(p)

    elif aln == "all":
        for p in [mafftdir, muscledir, clustaldir]:
            if not os.path.exists(p):
                os.mkdir(p)
        
    return mafftdir, muscledir, clustaldir, macsecdir, macseodir

def check_clean(indir, logpath):
    os.chdir(indir)
    
    prior_files = [f for f in os.listdir('.') if f.endswith(("clustalo_temp.fasta", "mafft_temp.fasta",
                                                                 "muscle_temp.fasta", "MAFFT_Aligned.fasta",
                                                                 "MUSCLE_Aligned.fasta", "CLUSTALO_Aligned.fasta",
                                                                 "_AA.fasta", "_NT.fasta"))]
    if prior_files:
        message = ("WARNING: Found {:,} files that appear to be from a previous run:".format(len(prior_files)))
        print("\n\n{}".format(message))
        write_log(logpath, message)
        for f in prior_files:
            print("\t{}".format(f))
            write_log(logpath, "\t{}".format(f))
        print("\n\tRemoving these files now.\n")
        for f in prior_files:
            os.remove(f)
            
def write_log(logpath, s):
    with open(logpath, 'a') as fh:
        fh.write("{}\n".format(s))
       
def main():
    tb = datetime.now()
    args = get_args()
    argd = vars(args)
    settings = ("Align.py settings:\n-i: {0}\n"
                    "-o: {1}\n"
                    "-a: {2}\n"
                    "--accurate: {3}\n"
                    "--threads: {4}\n"
                    "--mpath: {5}\n"
                    "--table: {6}\n"
                    "--mem: {7}\n"
                    "--pass_fail: {8}\n".format(argd["indir"], argd["outdir"],
                                                   argd["aln"], argd["accurate"],
                                                    argd["threads"], argd["mpath"],
                                                    argd["table"], argd["mem"],
                                                    argd["pass_fail"]))
        
    print("\n\n{}\n\n".format(settings))
    logpath = os.path.join(args.outdir, "Align.log")
    write_log(logpath, "Run executed: {}\n\n{}\n\n".format(datetime.now(), settings))
    
    mafftdir, muscledir, clustaldir, macsecdir, macseodir = make_dirs(args.outdir, args.aln, logpath)

    check_clean(args.indir, logpath)

    os.chdir(args.indir)
    flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
    if not flist:
        raise ValueError(("\n\n\nNo files with the extension .fa or .fasta were found in the "
                              "input directory:\n\t{}\n\n\n".format(args.indir)))
    
    commands = get_all_commands(flist, args.aln, args.accurate, args.threads, args.mpath, args.table, args.mem, args.pass_fail)
    
    for c in commands:
        print("\n\n{}".format(c))
        write_log(logpath, "{}: Executed: {}".format(datetime.now(), c))
        sp.call(c, shell=True)
        cleanup(args.aln, mafftdir, muscledir, clustaldir, macsecdir, macseodir, logpath)
    
    tf = datetime.now()
    print("\n\n{}".format("="*90))
    print("\nTotal time to create {2} aligments using {1}: {0} (H:M:S)\n".format(tf - tb, args.aln, len(commands)))
    print("{}\n\n".format("="*90))
        

if __name__ == '__main__':
    main()
