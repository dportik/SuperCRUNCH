''''
SuperCRUNCH: Align module

Usage: python Align.py  -i [directory with all fasta files] (REQUIRED)
                        -a [mafft, macse, muscle, clustalo] (REQUIRED)
                        --mpath [Full path to a macse.jar file] (Required for -a macse)
                        --mem [An integer for how much memory to assign to macse (in GB), default=1] (Optional for -a macse)
                        --table [Specifies translation table for macse (standard, vmtdna), default=standard] (Optional for -a macse)
                        --pass_fail (Optional for -a macse)
                        --accurate (Optional for -a mafft/macse)
                        --threads [integer] (Optional for mafft and clustalo)
						
    Align: Perform alignments for a directory of unaligned fasta files using mafft, macse, 
    muscle, or clustalo. 
    
    Empirical fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME. See below for
    special case with macse.

    Select analysis type with -a flag (options: mafft, macse, muscle, clustalo). If macse option is selected, user must
    provide the --mpath flag with the full path to a macse.jar file. Optional for macse are assigning
    a value for memory in GB (--mem, default 1GB) and to specify a different translation table
    (--table, default = standard code). The current options include: standard, vertmtdna, invertmtdna, yeastmtdna, plastid, 1-6, 9-16, 21-23. 
    The macse alignments can take a long time to complete
    so the elapsed time is shown after an alignment is produced for a given fasta file. Because
    macse will insert an ! at corrected bp locations, a cleaned fasta file (in which ! is replaced 
    by N) is also output.

    An additional feature of macse is to include a fasta file of reliable sequences
    (for example those that passed translation) and a fasta file of less reliable sequences
    that are suspected to contain errors, and align both simultaneously with different
    penalty parameters. To use this feature the --pass_fail flag can be used. However, to work
    the fasta files must follow this naming format:
    
    'prefix_Passed.fasta' - Reliable sequences fasta file.
    'prefix_Failed.fasta' - Unreliable sequences fasta file.
    
    The prefix portion of the name cannot contain any underscores, and should ideally just be
    the name of the gene/locus. 

    For other analyses, mafft, muscle, and clustalo must be installed in path (with identical 
    executable names as just written) for this script to work properly. These analyses will 
    run under the default settings (muscle) or the auto select settings (mafft,
    clustalo), but see below for an additional mafft, clustalo, and macse setting. The number of 
    threads can be specified for mafft and clustalo using the --threads flag.
    
    The optional flag --accurate can be used for mafft, clustalo, and macse v2.0+ to change the analysis 
    settings for increased accuracy. For mafft, this option changes the default from auto select
    to use FFT-NS-i ('mafft --retree 2 --maxiterate 1000') iterative search setting which may  
    result in higher accuracy (but slower run times). For clustalo, this de-selects the --auto
    option and enables --iter=5, in which the guide-tree and HMM each undergo 5 iterations, rather
    than only one. For macse, this adds the commands -local_realign_init 0.9
    -local_realign_dec 0.9 (both defaults = 0.5), which will slow down optimizations but increase
    alignment accuracy, sometimes considerably. If using macse v2.0+, these search settings will
    more closely resemble macse v1, for which both defaults = 1.0.
 
    Output files vary between aligners but will be moved to following output directories
    created in the main fasta directory:

    /Output_MACSE_Alignments:
           [fasta name]_AA.fasta - The amino acid alignment.
           [fasta name]_NT.fasta - The nucleotide alignment from the translation alignment. 
           [fasta name]_NT_Cleaned.fasta - Same as above, but the ! characters inserted by macse
                                           are replaced with Ns instead.
                                           
    /Output_MAFFT_Alignments:
           [fasta name]_MAFFT_Aligned.fasta - The nucleotide alignment produced by mafft.
           
    /Output_MUSCLE_Alignments:
           [fasta name]_MUSCLE_Aligned.fasta - The nucleotide alignment produced by muscle.
           
    /Output_CLUSTALO_Alignments:
           [fasta name]_CLUSTALO_Aligned.fasta - The nucleotide alignment produced by clustalo.
    
-------------------------
For Python 2.7
Python modules required:
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
January 2019
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
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Align: Perform alignments for a directory of unaligned fasta files using mafft, macse, 
    muscle, or clustalo. 
    
    Empirical fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME. See below for
    special case with macse.

    Select analysis type with -a flag (options: mafft, macse, muscle, clustalo). If macse option is selected, user must
    provide the --mpath flag with the full path to a macse.jar file. Optional for macse are assigning
    a value for memory in GB (--mem, default 1GB) and to specify a different translation table
    (--table, default = standard code). The current options include: standard, vertmtdna, invertmtdna, yeastmtdna, plastid, 1-6, 9-16, 21-23. 
    The macse alignments can take a long time to complete
    so the elapsed time is shown after an alignment is produced for a given fasta file. Because
    macse will insert an ! at corrected bp locations, a cleaned fasta file (in which ! is replaced 
    by N) is also output.

    An additional feature of macse is to include a fasta file of reliable sequences
    (for example those that passed translation) and a fasta file of less reliable sequences
    that are suspected to contain errors, and align both simultaneously with different
    penalty parameters. To use this feature the --pass_fail flag can be used. However, to work
    the fasta files must follow this naming format:
    
    'prefix_Passed.fasta' - Reliable sequences fasta file.
    'prefix_Failed.fasta' - Unreliable sequences fasta file.
    
    The prefix portion of the name cannot contain any underscores, and should ideally just be
    the name of the gene/locus. 

    For other analyses, mafft, muscle, and clustalo must be installed in path (with identical 
    executable names as just written) for this script to work properly. These analyses will 
    run under the default settings (muscle) or the auto select settings (mafft,
    clustalo), but see below for an additional mafft, clustalo, and macse setting. The number of 
    threads can be specified for mafft and clustalo using the --threads flag.
    
    The optional flag --accurate can be used for mafft, clustalo, and macse v2.0+ to change the analysis 
    settings for increased accuracy. For mafft, this option changes the default from auto select
    to use FFT-NS-i ('mafft --retree 2 --maxiterate 1000') iterative search setting which may  
    result in higher accuracy (but slower run times). For clustalo, this de-selects the --auto
    option and enables --iter=5, in which the guide-tree and HMM each undergo 5 iterations, rather
    than only one. For macse, this adds the commands -local_realign_init 0.9
    -local_realign_dec 0.9 (both defaults = 0.5), which will slow down optimizations but increase
    alignment accuracy, sometimes considerably. If using macse v2.0+, these search settings will
    more closely resemble macse v1, for which both defaults = 1.0.

    Output files vary between aligners but will be moved to output directories
    created in the main fasta directory.
    
    DEPENDENCIES: Python: BioPython; Executables in path: mafft, muscle, clustalo (with
    names identical to those listed here); macse jar file required.
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files. Follow labeling format: NAME.fasta")
    parser.add_argument("-a", "--aln", required=True, choices=["mafft","macse","muscle", "clustalo","all"], help="REQUIRED: Specify whether alignment is by mafft, macse, muscle, or clustalo. If macse must provide flags --mpath with full path to macse jar file and --table with translation table option. Selecting 'all' will run mafft, muscle, and clustalo sequentially.")
    parser.add_argument("--mpath", default=None, help="MACSE REQUIRED: Full path to a macse.jar file.")
    parser.add_argument("--table", default="standard", choices=["standard","vertmtdna","invertmtdna","yeastmtdna","plastid","1","2","3","4","5","6","9","10","11","12","13","14","15","16","21","22","23"], help="MACSE REQUIRED: Specifies translation table for macse.")
    parser.add_argument("--mem", default=None, help="MACSE OPTIONAL: An integer for how much memory to assign to macse (in GB), default=1.")
    parser.add_argument("--pass_fail", action='store_true', help="MACSE OPTIONAL: Specifies macse to use two fasta files (one with seqs passing translation, and one with those that failed) for dual file alignment. Follow labeling format: NAME_Passed.fasta, NAME_Failed.fasta")
    parser.add_argument("--accurate", action='store_true', help="MACSE/MAFFT OPTIONAL: Specifies mafft, clustalo, or macse to use more thorough search settings.")
    parser.add_argument("--threads", default=None, help="OPTIONAL: Specifies number of threads to use in mafft or clustalo.")
    return parser.parse_args()

def table_dict(symbol):
    '''
    Get MACSE specific translation table value from table choice.
    '''
    table_dictv1 = {"standard":"1","vertmtdna":"2","invertmtdna":"5","yeastmtdna":"3","plastid":"11","1":"1","2":"2","3":"3","4":"4","5":"5","6":"6","9":"9","10":"10","11":"11","12":"12","13":"13","14":"14","15":"15","16":"16","21":"21","22":"22","23":"23"}
    table_dictv2 = {"standard":"01_The_Standard_Code", "vertmtdna":"02_The_Vertebrate_Mitochondrial_Code", "invertmtdna":"05_The_Invertebrate_Mitochondrial_Code", "yeastmtdna":"03_The_Yeast_Mitochondrial_Code", "plastid":"11_The_Bacterial_Archaeal_and_Plant_Plastid_Code", "1":"01_The_Standard_Code", "2":"02_The_Vertebrate_Mitochondrial_Code", "3":"03_The_Yeast_Mitochondrial_Code", "4":"04_The_Mold_Protozoan_and_Coelenterate_Mitochondrial_Code_and_the_Mycoplasma_Spiroplasma_Code", "5":"05_The_Invertebrate_Mitochondrial_Code", "6":"06_The_Ciliate_Dasycladacean_and_Hexamita_Nuclear_Code", "9":"09_The_Echinoderm_and_Flatworm_Mitochondrial_Code", "10":"10_The_Euplotid_Nuclear_Code", "11":"11_The_Bacterial_Archaeal_and_Plant_Plastid_Code", "12":"12_The_Alternative_Yeast_Nuclear_Code", "13":"13_The_Ascidian_Mitochondrial_Code", "14":"14_The_Alternative_Flatworm_Mitochondrial_Code", "15":"15_Blepharisma_Nuclear_Code", "16":"16_Chlorophycean_Mitochondrial_Code", "21":"21_Trematode_Mitochondrial_Code", "22":"22_Scenedesmus_obliquus_mitochondrial_Code", "23":"23_Thraustochytrium_Mitochondrial_Code"}
    table_val = table_dictv1[symbol]
    return table_val

def directory_macse_aln(in_dir, mpath, mem, incode, acc, flist):
    '''
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the macse_align function
    for each file found. Moves all output files to the output
    directory: /Output_MACSE_Alignments
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tBeginning MACSE alignments"
    print "--------------------------------------------------------------------------------------"
    
    tcode = table_dict(incode)
    
    os.chdir(in_dir)
    out_dir = "Output_MACSE_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    for fasta in flist:
        macse_align(fasta, mpath, mem, tcode, acc)
        
    output = [f for f in os.listdir('.') if f.endswith(('_AA.fasta', '_NT.fasta'))]
    for o in output:
        shutil.move(o, out_dir)
                    
    #clean NT.fasta files (replace !s with Ns)
    os.chdir(out_dir)
    out_list = [f for f in os.listdir('.') if f.endswith("_NT.fasta")]
    for f in out_list:
        cleanlabel = "{}_NT_Cleaned.fasta".format((f.split('_NT')[0]))
        with open(cleanlabel, 'a') as fh_out:
            with open(f, 'r') as fh_in:
                for line in fh_in:
                    line = line.replace('!','N')
                    fh_out.write(line)
                                        
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tFinished MACSE alignments"
    print "--------------------------------------------------------------------------------------\n\n"

def macse_align(fasta_file, mpath, mem, tcode, acc):
    print "\n\nPerforming MACSE alignment for {}\n\n".format(fasta_file)
    t_begin = datetime.now()
    #gcdef codes: (1:Standard; 2:Vert mtDNA)
    #macse can automatically name output files here so no need to define prefix
    if mem is not None:
        if acc is True:
            call_string = "java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -local_realign_init 0.9 -local_realign_dec 0.9 ".format(mem, mpath, tcode, fasta_file)
        else:
            call_string = "java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} ".format(mem, mpath, tcode, fasta_file)
    else:
        if acc is True:
            call_string = "java -jar -Xmx1g {0} -prog alignSequences -gc_def {1} -seq {2} -local_realign_init 0.9 -local_realign_dec 0.9 ".format(mpath, tcode, fasta_file)
        else:
            call_string = "java -jar -Xmx1g {0} -prog alignSequences -gc_def {1} -seq {2} ".format(mpath, tcode, fasta_file)
        
    print call_string
    proc = sp.call(call_string, shell=True)
    t_finish = datetime.now()
    elapsed = t_finish - t_begin
    print "Total alignment time: {0} (H:M:S)\n\n".format(elapsed)


def split_name(string, index, delimiter):
    index = int(index)
    name = string.split(delimiter)[index]
    return name

def pass_fail_finder(in_dir):
    prefix_set = set()
    os.chdir(in_dir)
    for f in os.listdir('.'):
    	if f.endswith('_Passed.fasta'):
        	name = split_name(f, 0, '_Passed')
        	prefix_set.add(name)
    prefix_list = list(prefix_set)
    prefix_list.sort()
    return prefix_list
    
def directory_macse_aln_pass_fail(in_dir, mpath, mem, incode, acc):
    '''
    Iterates over files in a directory to locate those with
    same prefix but different extensions '_Passed.fasta' and
    '_Failed.fasta' and executes the macse_align function using
    the -seq and -seqlr options for each file pair found.
    Moves all output files to the output directory:
    /Output_MACSE_Alignments
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tBeginning paired MACSE alignments"
    print "--------------------------------------------------------------------------------------"

    tcode = table_dict(incode)
    
    prefix_list = pass_fail_finder(in_dir)

    os.chdir(in_dir)
    
    out_dir = "Output_MACSE_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for p in prefix_list:
        f1 = "{}_Passed.fasta".format(p)
        f2 = "{}_Failed.fasta".format(p)
        macse_align_pass_fail(f1, f2, mpath, mem, tcode, acc)
        
    output = [f for f in os.listdir('.') if f.endswith(("_AA.fasta", "_NT.fasta"))]
    for o in output:
        shutil.move(o, out_dir)
                    
    #clean NT.fasta files (replace !s with Ns)
    os.chdir(out_dir)
    out_list = [f for f in os.listdir('.') if f.endswith("_NT.fasta")]
    for f in out_list:
        cleanlabel = "{}_NT_Cleaned.fasta".format((f.split('_NT')[0]))
        with open(cleanlabel, 'a') as fh_out:
            with open(f, 'r') as fh_in:
                for line in fh_in:
                    line = line.replace('!','N')
                    fh_out.write(line)
                    
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tFinished paired MACSE alignments"
    print "--------------------------------------------------------------------------------------\n\n"
    
def macse_align_pass_fail(f1, f2, mpath, mem, tcode, acc):
    print "\n\nPerforming MACSE alignment for {0} and {1}\n\n".format(f1,f2)
    t_begin = datetime.now()
    
    #macse automatically names output files here so no need to define prefix
    if mem is not None:
        if acc is True:
            call_string = "java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -seq_lr {4} -local_realign_init 0.9 -local_realign_dec 0.9 ".format(mem, mpath, tcode, f1, f2)
        else:
            call_string = "java -jar -Xmx{0}g {1} -prog alignSequences -gc_def {2} -seq {3} -seq_lr {4} ".format(mem, mpath, tcode, f1, f2)
    else:
        if acc is True:
            call_string = "java -jar -Xmx1g {0} -prog alignSequences -gc_def {1} -seq {2} -seq_lr {3} -local_realign_init 0.9 -local_realign_dec 0.9 ".format(mpath, tcode, f1, f2)
        else:
            call_string = "java -jar -Xmx1g {0} -prog alignSequences -gc_def {1} -seq {2} -seq_lr {3} ".format(mpath, tcode, f1, f2)
        
    print call_string
    proc = sp.call(call_string, shell=True)
    t_finish = datetime.now()
    elapsed = t_finish - t_begin
    print "Total alignment time: {0} (H:M:S)\n\n".format(elapsed)
    
def directory_mafft_aln(in_dir, acc, threads, flist):
    '''
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the mafft_align function
    for each file found. Moves all output files to the output
    directory: /Output_MAFFT_Alignments
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tBeginning MAFFT alignments"
    print "--------------------------------------------------------------------------------------"
    os.chdir(in_dir)

    out_dir = "Output_MAFFT_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    for fasta in flist:
        mafft_align(fasta, acc, threads)
        
    output = [f for f in os.listdir('.') if f.endswith("MAFFT_Aligned.fasta")]
    for o in output:
        try:
            shutil.move(o, out_dir)
        except:
            print "{} was not moved. It is already in the MAFFT alignment output directory, you may want to check this.".format(o)
                    
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tFinished MAFFT alignments"
    print "--------------------------------------------------------------------------------------\n\n"

def mafft_align(fasta_file, acc, threads):
    print "\n\nPerforming MAFFT alignment for {}\n\n".format(fasta_file)
    names = fasta_file.split('.')
    
    if threads is None:
        threads = 1
    if acc is False:
    	call_string = "mafft --thread {2} --auto {0} > {1}_mafft_temp.fasta".format(fasta_file, names[0], threads)
    else:
    	call_string = "mafft --thread {2} --retree 2 --maxiterate 1000 {0} > {1}_mafft_temp.fasta".format(fasta_file, names[0], threads)
        
    print call_string, '\n'
    proc = sp.call(call_string, shell=True)

    fasta_dict = SeqIO.index("{0}_mafft_temp.fasta".format(names[0]), "fasta")
    out_fasta = "{0}_MAFFT_Aligned.fasta".format(names[0])
    with open(out_fasta, 'a') as fh_out_fasta:
        for record in fasta_dict:
            newseq = fasta_dict[record].seq.upper()
            fh_out_fasta.write( ">{}\n{}\n".format(fasta_dict[record].description, newseq))
    os.remove("{0}_mafft_temp.fasta".format(names[0]))
    
def directory_muscle_aln(in_dir, flist):
    '''
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the mafft_align function
    for each file found. Moves all output files to the output
    directory: /Output_MUSCLE_Alignments
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tBeginning MUSCLE alignments"
    print "--------------------------------------------------------------------------------------"
    os.chdir(in_dir)

    out_dir = "Output_MUSCLE_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    for fasta in flist:
        muscle_align(fasta)
        
    output = [f for f in os.listdir('.') if f.endswith("muscle_temp.fasta")]
    for o in output:
        shutil.move(o, out_dir)
                    
    os.chdir(out_dir)
    out_list = [f for f in os.listdir('.') if f.endswith("muscle_temp.fasta")]
    for f in out_list:
        fasta_dict = SeqIO.index(f, "fasta")
        out_fasta = "{0}_MUSCLE_Aligned.fasta".format(f.split('_muscle_temp')[0])
        with open(out_fasta, 'a') as fh_out_fasta:
            for record in fasta_dict:
                newseq = fasta_dict[record].seq.upper()
                fh_out_fasta.write( ">{}\n{}\n".format(fasta_dict[record].description, newseq))
    for f in out_list:
        os.remove(f)
            
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tFinished MUSCLE alignments"
    print "--------------------------------------------------------------------------------------\n\n"

def muscle_align(fasta_file):
    print "\n\nPerforming MUSCLE alignment for {}\n\n".format(fasta_file)
    names = fasta_file.split('.')
    
    call_string = "muscle -in {0} -out {1}_muscle_temp.fasta".format(fasta_file, names[0])
    
    print call_string, '\n'
    proc = sp.call(call_string, shell=True)

def directory_clustalo_aln(in_dir, acc, threads, flist):
    '''
    Iterates over files in a directory to locate those with
    extension '.fasta' and executes the clustalo_align function
    for each file found. Moves all output files to the output
    directory: /Output_CLUSTALO_Alignments
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tBeginning Clustal-Omega alignments"
    print "--------------------------------------------------------------------------------------"
    os.chdir(in_dir)

    out_dir = "Output_CLUSTALO_Alignments"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    for fasta in flist:
        clustalo_align(fasta, acc, threads)
        
    output = [f for f in os.listdir('.') if f.endswith("clustalo_temp.fasta")]
    for o in output:
        shutil.move(o, out_dir)
                    
    os.chdir(out_dir)
    out_list = [f for f in os.listdir('.') if f.endswith("clustalo_temp.fasta")]
    for f in out_list:
        fasta_dict = SeqIO.index(f, "fasta")
        out_fasta = "{0}_CLUSTALO_Aligned.fasta".format(f.split('_clustalo_temp')[0])
        with open(out_fasta, 'a') as fh_out_fasta:
            for record in fasta_dict:
                newseq = fasta_dict[record].seq.upper()
                fh_out_fasta.write( ">{}\n{}\n".format(fasta_dict[record].description, newseq))
    for f in out_list:
        os.remove(f)
            
    print "\n\n--------------------------------------------------------------------------------------"
    print "\t\t\tFinished Clustal-Omega alignments"
    print "--------------------------------------------------------------------------------------\n\n"

def clustalo_align(fasta_file, acc, threads):
    print "\n\nPerforming Clustal-Omega alignment for {}\n\n".format(fasta_file)
    names = fasta_file.split('.')
    
    if threads is None:
        threads = 1
    if acc is False:
        call_string = "clustalo -i {0} -o {1}_clustalo_temp.fasta --auto -v --threads={2} --output-order=tree-order --force".format(fasta_file, names[0], threads)
    else:
        call_string = "clustalo -i {0} -o {1}_clustalo_temp.fasta --full --full-iter --iter=5 -v --threads={2} --cluster-size=500 --output-order=tree-order --force".format(fasta_file, names[0], threads)
        
    print call_string, '\n'
    proc = sp.call(call_string, shell=True)

#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    
    os.chdir(args.in_dir)
    finitial = sorted([f for f in os.listdir('.') if f.endswith((".fasta", ".fa"))])
    flist = sorted([f for f in finitial if not f.endswith(("clustalo_temp.fasta","mafft_temp.fasta","muscle_temp.fasta", "MAFFT_Aligned.fasta", "MUSCLE_Aligned.fasta", "CLUSTALO_Aligned.fasta"))])
    
    if args.aln == "mafft":
        tb = datetime.now()
        directory_mafft_aln(args.in_dir, args.accurate, args.threads, flist)
        tf = datetime.now()
        te = tf - tb
        print "Total time for all alignments using {0}: {1} (H:M:S)\n\n".format(args.aln,te)
        
    elif args.aln == "macse":
        if args.pass_fail is True:
            tb = datetime.now()
            directory_macse_aln_pass_fail(args.in_dir, args.mpath, args.mem, args.table, args.accurate)
            tf = datetime.now()
            te = tf - tb
            print "Total time for all paired translation alignments using {0}: {1} (H:M:S)\n\n".format(args.aln,te)
            
        elif args.pass_fail is False:
            tb = datetime.now()
            directory_macse_aln(args.in_dir, args.mpath, args.mem, args.table, args.accurate, flist)
            tf = datetime.now()
            te = tf - tb
            print "Total time for all translation alignments using {0}: {1} (H:M:S)\n\n".format(args.aln,te)
            
    elif args.aln == "muscle":
        tb = datetime.now()
        directory_muscle_aln(args.in_dir, flist)
        tf = datetime.now()
        te = tf - tb
        print "Total time for all alignments using {0}: {1} (H:M:S)\n\n".format(args.aln,te)
        
    elif args.aln == "clustalo":
        tb = datetime.now()
        directory_clustalo_aln(args.in_dir, args.accurate, args.threads, flist)
        tf = datetime.now()
        te = tf - tb
        print "Total time for all alignments using {0}: {1} (H:M:S)\n\n".format(args.aln,te)
        
    elif args.aln == "all":
        tb = datetime.now()

        tbmaf = datetime.now()
        directory_mafft_aln(args.in_dir, args.accurate, args.threads, flist)
        tfmaf = datetime.now()
        temaf = tfmaf - tbmaf

        tbclu = datetime.now()
        directory_clustalo_aln(args.in_dir, args.accurate, args.threads, flist)
        tfclu = datetime.now()
        teclu = tfclu - tbclu
        
        tbmus = datetime.now()
        directory_muscle_aln(args.in_dir, flist)
        tfmus = datetime.now()
        temus = tfmus - tbmus
        
        tf = datetime.now()
        te = tf - tb
        
        print "Total time for all alignments using MAFFT: {0} (H:M:S)\n\n".format(temaf)
        print "Total time for all alignments using CLUSTAL-O: {0} (H:M:S)\n\n".format(teclu)
        print "Total time for all alignments using MUSCLE: {0} (H:M:S)\n\n".format(temus)      
        print "Total time to run all aligners: {0} (H:M:S)\n\n".format(te)
        
    else:
        print "Something went wrong!"

if __name__ == '__main__':
    main()
