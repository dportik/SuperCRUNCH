''''
SuperCRUNCH: Coding_Translation_Tests module

    Coding_Translation_Tests: This function processes an unaligned fasta file. 
    If necessary, it will adjust sequences with 'N's so they are divisible by 
    three to complete the final 'codon'. Sequences are translated in all forward 
    frames to check for the presence of stop codons. If stop codons are detected
    in all frames, it will next check to see if there is only one stop codon and 
    it is present in the final two codon positions of the sequence. If the 
    --rc flag is included the translation will also be performed for the reverse 
    complement, however if your sequences are all correctly oriented this is not 
    recommended. The final sequence is written such that the first base represents 
    the first codon position (based on the detected reading frame) and the sequence 
    is divisible by three (to ensure a complete final codon position).

    The translation table should be specified with the --table flag. All NCBI 
    translation table options are available, and can be selected using integers 
    or the shortcut terms provided. If the --table flag is omitted, the default 
    will be to use the Standard translation table. 

    Sequences that have more than one stop codon in all frames fail the translation 
    test. Since no reading frame is better than any other, the sequence will be 
    written as is (adjusted for length to be divisible by three) to the relevant 
    files, noted by an '*' appended at the end of the sequence description line.

    If your dataset contains a mix of loci that need to be translated 
    differently (e.g., coding mtDNA vs. coding nucDNA), the 
    --onlyinclude flag can be used. This requires the full path to a text file 
    containing the names of the fasta files (one file name per line) to be 
    processed for a particular run. This allows you to run this module with 
    particular settings for a subset of the files in the input directory.

    Perhaps the most important aspects of this script are the adjusting of sequences 
    to the first codon position and the completion of the final codon. This is 
    required for some alignment programs.
    
    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

    All output files are written to relevant directories in the specified
    output directory.

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
import os
import subprocess as sp
import shutil
import argparse
import re
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Coding_Translation_Tests: This function processes an unaligned fasta file. 
    If necessary, it will adjust sequences with 'N's so they are divisible by 
    three to complete the final 'codon'. Sequences are translated in all forward 
    frames to check for the presence of stop codons. If stop codons are detected
    in all frames, it will next check to see if there is only one stop codon and 
    it is present in the final two codon positions of the sequence. If the 
    --rc flag is included the translation will also be performed for the reverse 
    complement, however if your sequences are all correctly oriented this is not 
    recommended. The final sequence is written such that the first base represents 
    the first codon position (based on the detected reading frame) and the sequence 
    is divisible by three (to ensure a complete final codon position).

    The translation table should be specified with the --table flag. All NCBI 
    translation table options are available, and can be selected using integers 
    or the shortcut terms provided. If the --table flag is omitted, the default 
    will be to use the Standard translation table. 

    Sequences that have more than one stop codon in all frames fail the translation 
    test. Since no reading frame is better than any other, the sequence will be 
    written as is (adjusted for length to be divisible by three) to the relevant 
    files, noted by an '*' appended at the end of the sequence description line.

    If your dataset contains a mix of loci that need to be translated 
    differently (e.g., coding mtDNA vs. coding nucDNA), the 
    --onlyinclude flag can be used. This requires the full path to a text file 
    containing the names of the fasta files (one file name per line) to be 
    processed for a particular run. This allows you to run this module with 
    particular settings for a subset of the files in the input directory.

    Perhaps the most important aspects of this script are the adjusting of sequences 
    to the first codon position and the completion of the final codon. This is 
    required for some alignment programs.
    
    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

    All output files are written to relevant directories in the specified
    output directory.
    
    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to "
                            "write output files.")
    
    parser.add_argument("--table",
                            default="standard",
                            choices=["standard", "vertmtdna", "invertmtdna", "yeastmtdna",
                                         "plastid", "1", "2", "3", "4", "5", "6", "7", "8",
                                         "9", "10", "11", "12", "13", "14", "15", "16", "17",
                                         "18", "19", "20", "21", "22", "23", "24", "25", "26", 
                                         "27", "28", "29", "30", "31"],
                            help="REQUIRED: Specifies translation table to use for all files.")
    
    parser.add_argument("--rc",
                            action='store_true',
                            help="OPTIONAL: In addition to all forward reading frames, examine "
                            "all reverse reading frames during translation tests.")

    parser.add_argument("--onlyinclude",
                            required=False,
                            help="OPTIONAL: The full path to a text file containing the "
                            "names of the fasta files to process for this run.")
    
    parser.add_argument("--quiet",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Show less output while running (useful when filtering "
                            "many loci).")
   
    return parser.parse_args()

def translation_val(tcode):
    """
    Convert shortcut translation table terms into the true NCBI
    nomenclature using a dictionary structure. Term needed for 
    translation with biopython.
    """
    tdict = {"standard":"Standard", "vertmtdna":"Vertebrate Mitochondrial",
                 "invertmtdna":"Invertebrate Mitochondrial", "yeastmtdna":"Yeast Mitochondrial",
                 "plastid":"11", "1":1, "2":2, "3":3, "4":4, "5":5, "6":6, "7":7, "8":8,
                 "9":9, "10":10, "11":11, "12":12, "13":13, "14":14, "15":15, "16":16,
                 "17":17, "18":18, "19":19, "20":20, "21":21, "22":22, "23":23, "24":24,
                 "25":25, "26":26, "27":27, "28":28, "29":29, "30":30, "31":31}
        
    tval = tdict[tcode]
    
    return tval

def orf_adjust(f, tcode, rc, quiet):
    """
    Function to adjust all sequences in file (f) to the proper
    reading frame (based on tcode translation), if possible.

    Three output files are produced:
    [input fasta name]_All.fasta - Contains all sequences, pass and fail.
    [input fasta name]_Passed.fasta - Contains only sequences that passed translation.
    [input fasta name]_Failed.fasta - Contains only sequences that failed translation.
    """

    #initiate counts of sequences passing or failing translations
    tpass, tfail = int(0), int(0)
    
    #labels for output files
    prefix = f.split('/')[-1].split('.')[0]
    fasta_all = "{}_All.fasta".format(prefix)
    fasta_pass = "{}_Passed.fasta".format(prefix)
    fasta_fail = "{}_Failed.fasta".format(prefix)

    #index fasta file and iterate
    fasta_dict = SeqIO.index(f, "fasta")

    #initiate loop to check every sequence
    for record in fasta_dict:
        
        #adjust sequence length to be divisible by three for translation
        #check, otherwise biopython gets angry if final 'codon' incomplete
        #length adjusted by adding N's to end
        newseq = fasta_dict[record].seq
        n = Seq("N")
        if len(newseq) % 3 == 0:
            pass
        
        elif len(newseq) %3 != 0:
            newseq = newseq + n
            if len(newseq) % 3 == 0:
                pass
            
            elif len(newseq) % 3 != 0:
                newseq = newseq + n
        
        #perform translation with table selected by user
        #translate to all forward frames
        #for each frame considered, adjust length to a multiple of 3
        frame1 = newseq[:].translate(table=tcode)

        #add 1 N, slice sequence to exclude first base position
        f2 = newseq+n
        frame2 = f2[1:].translate(table=tcode)

        #add another N, slice sequence to exclude first+second base positions
        f3 = f2+n
        frame3 = f3[2:].translate(table=tcode)
        
        #reverse complement and translate to all frames
        #as above, adjust to keep multiples of 3
        rcseq = newseq.reverse_complement()
        rcframe1 = rcseq[:].translate(table=tcode)
        
        rf2 = rcseq+n
        rcframe2 = rf2[1:].translate(table=tcode)
        
        rf3 = rf2+n
        rcframe3 = rf3[2:].translate(table=tcode)

        #count stop codons present in each frame
        f1stops = frame1.count('*')
        f2stops = frame2.count('*')
        f3stops = frame3.count('*')
        rcf1stops = rcframe1.count('*')
        rcf2stops = rcframe1.count('*')
        rcf3stops = rcframe1.count('*')

        #create sublists with relevant information to iterate over
        #for each frame - [stop codon count, sequence object in that
        #frame, frame string label, seq slice divisible by 3 based on frame]
        
        #if rc option is included, use the rc as well
        if rc is False:
            frames_list = [[f1stops,frame1,"Frame 1",newseq[:]],
                      [f2stops,frame2,"Frame 2",f2[1:]],
                      [f3stops,frame3,"Frame 3",f3[2:]]]
                
        elif rc is True:
            frames_list = [[f1stops,frame1,"Frame 1",newseq[:]],
                      [f2stops,frame2,"Frame 2",f2[1:]],
                      [f3stops,frame3,"Frame 3",f3[2:]],
                      [rcf1stops,rcframe1,"RC Frame 1",rcseq[:]],
                      [rcf2stops,rcframe2,"RC Frame 2",rf2[1:]],
                      [rcf3stops,rcframe3,"RC Frame 3",rf3[2:]]]

        #Decision tree for sequence selection
        #first: we search all frames to see if one has 0 stop codons, if so
        #       we use that sequence frame
        #second: we search all frames to see if a single stop codon occurs
        #        in the final two positions, if so we use that frame
        #third: if above two searches fail we just use length-adjusted
        #       sequence in whatever frame it came in

        #set match to False to begin search
        match = False

        #this will check for 0 stop codons in frames
        for sublist in frames_list:
            if match is False:        
                #check if frame has no stop codons
                if sublist[0] == int(0):
                    tpass+=1
                    if not quiet:
                        print("\t{0}: Translation {1}".format(fasta_dict[record].id, sublist[2]))
                    
                    #make sure that any final codon composed entirely of Ns is removed from seq
                    #use regular expression to find and remove
                    write_seq = re.sub('NNN$', '', str(sublist[3]))

                    #write to output files
                    with open(fasta_all, 'a') as fh_all:
                        fh_all.write(">{0}\n{1}\n".format(fasta_dict[record].description, write_seq) )
                        
                    with open(fasta_pass, 'a') as fh_pass:
                        fh_pass.write(">{0}\n{1}\n".format(fasta_dict[record].description, write_seq) )

                    #set match to true to skip next steps in decision tree
                    match = True

        #this will check for a stop codon in the final two positions
        if match is False:
            for sublist in frames_list:
                #tricky part - somewhere in this search we may hit a
                #good seq frame and turn match to true, so check
                #every time we move to another seq frame and break
                #when appropriate
                if match is False:
                    #check if frame has one stop codon and it is in the last two positions of the sequence
                    if sublist[0] == int(1) and '*' in sublist[1][-3:]:
                        tpass+=1
                        if not quiet:
                            print("\t{0}: Translation {1} with single tailing stop codon"
                                      .format(fasta_dict[record].id, sublist[2]))
                            for i in frames_list:
                                print("\t\t\t{} stop codons: {}".format(i[2], i[0]))
                            
                        #remove last two codons (including the detected stop codon)
                        #**turns out these tailing N's can interfere with alignment and other
                        #   downstream steps because entire columns can contain N's
                        write_seq = str(sublist[3])[:-6]

                        #write to output files
                        with open(fasta_all, 'a') as fh_all:
                            fh_all.write( ">{0}\n{1}\n".format(fasta_dict[record].description, write_seq))
                            
                        with open(fasta_pass, 'a') as fh_pass:
                            fh_pass.write( ">{0}\n{1}\n".format(fasta_dict[record].description, write_seq))
                            
                        #set match to true to break search
                        match = True

        #if all else fails...
        if match is False:
            if not quiet:
                print("\t{}: WARNING - Failed translation test".format(fasta_dict[record].id))
                for i in frames_list:
                    print("\t\t\t{0} stop codons: {1}".format(i[2], i[0]))
                    print('\t\t\tAdding * to end of the sequence description for {0}.'
                              .format(fasta_dict[record].id))
            
            tfail+=1
            
            #make sure that any final codon composed entirely of Ns is removed from seq
            write_seq = re.sub('NNN$', '', str(newseq))
            
            with open(fasta_all, 'a') as fh_all:
                fh_all.write(">{0}\n{1}\n".format("{} *".format(fasta_dict[record].description), write_seq))
                
            with open(fasta_fail, 'a') as fh_fail:
                fh_fail.write(">{0}\n{1}\n".format("{} *".format(fasta_dict[record].description), write_seq))

    summary = [prefix, tpass, tfail]
    
    return summary
                        
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
    pdir = os.path.join(curpath, "Translation-Passed-Seqs")
    if not os.path.exists(pdir):
        os.mkdir(pdir)
        
    fdir = os.path.join(curpath, "Translation-Failed-Seqs")
    if not os.path.exists(fdir):
        os.mkdir(fdir)
        
    adir = os.path.join(curpath, "Translation-All-Seqs")
    if not os.path.exists(adir):
        os.mkdir(adir)

    return pdir, fdir, adir

def write_log(results, outdir):
    os.chdir(outdir)
    
    outlog = "Log_Sequences_Filtered.txt"
    with open(outlog, 'a') as fh:
        fh.write("Locus\tSeqs_Passed\tSeqs_Failed\n")
        for r in results:
            fh.write("{}\t{}\t{}\n".format(r[0], r[1], r[2]))

def cleanup(pdir, fdir, adir):
    """
    Moves relevant output files to their output directories. 
    Deletes temporary or empty files.
    """
    [os.remove(f) for f in os.listdir('.') if os.stat(f).st_size == 0 and f.endswith('.fasta')]
    [shutil.move(f, pdir) for f in os.listdir('.') if f.endswith('Passed.fasta')]
    [shutil.move(f, fdir) for f in os.listdir('.') if f.endswith('Failed.fasta')]
    [shutil.move(f, adir) for f in os.listdir('.') if f.endswith('All.fasta')]
    
def main():
    args = get_args()
    tb = datetime.now()
    
    if args.outdir == '.':
        raise ValueError('\n\n***Please provide full path to output directory!\n')
    if args.indir == '.':
        raise ValueError('\n\n***Please provide full path to input directory!\n')
    
    pdir, fdir, adir = make_dirs(args.outdir)
            
    os.chdir(args.indir)
    
    if not args.onlyinclude:
        flist = sorted([os.path.abspath(f) for f in os.listdir('.')
                            if f.endswith(('.fasta', '.fa'))])
        print("\nFound {:,} fasta files to perform length adjustments and translation tests.".format(len(flist)))
        
    elif args.onlyinclude:
        with open(args.onlyinclude, 'r') as fh:
            incl = [l.strip() for l in fh]
        flist = sorted([os.path.abspath(f) for f in os.listdir('.')
                            if f.endswith(('.fasta', '.fa')) and f in incl])
        print("\nFound {:,} fasta files to perform length adjustments and translation tests, based on list of files to include.".format(len(flist)))
        
    if args.table is None:
        tcode = "Standard"
    else:
        tcode = translation_val(args.table)

    print("\n\nPerforming all translations using table: {}\n\n".format(tcode))

    results = []
    for f in flist:
        b = datetime.now()
        if not args.quiet:
            print("\n--------------------------------------------------------------------------------------\n")
        print("Adjusting reading frames of sequences in {}".format(f))
        fresults = orf_adjust(f, tcode, args.rc, args.quiet)
        results.append(fresults)
        cleanup(pdir, fdir, adir)
        f = datetime.now()
        e = f - b
        print("\tFinished. Elapsed time: {0} (H:M:S)\n".format(e))
        
    write_log(results, args.outdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal time to translate all files: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    

if __name__ == '__main__':
    main()
