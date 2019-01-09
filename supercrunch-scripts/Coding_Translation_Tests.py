''''
SuperCRUNCH: Coding_Translation_Tests module

Usage: python Coding_Translation_Tests.py -i [directory with all fasta files] (REQUIRED)
                                      --table ["standard","vertmtdna","invertmtdna","yeastmtdna","plastid","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"] (REQUIRED)
                                      --rc (Optional; also translate sequences in reverse complement)

    Coding_Translation_Tests: This function processes an unaligned fasta file. If necessary, it will adjust sequences with 'N's
    so they are divisible by 3 to complete the final codon. Sequences are translated in all
    forward frames to check for the presence of stop codons. If stop codons are detected
    in all frames, it will next check to see if there is only one stop codon and it is present in the
    final two codon positions of the sequence. If the --rc flag is included the translation will also
    be performed for the reverse complement, however if your sequences are all correctly oriented this
    is not recommended. The final sequence is written such that the first base represents the 
    first codon position (based on the detected reading frame) and the sequence is divisible 
    by three to ensure complete codons.

    The translation table should be specified with the --table flag. All NCBI translation table options 
    are available, and can be selected using integers or the shortcut terms provided. If the --table 
    flag is omitted, the default will be to use the Standard translation table. 

    Sequences that have more than one stop codon in all frames fail the translation test. Since no
    reading frame is better than any other, the sequence will be written as is (adjusted for length
    to be divisible by three) to the relevant files, noted by an '*' appended at the end of the sequence
    description line.

    The most important aspects of this script are the adjusting of sequences to the first codon position
    and the completion of the final codon. This is required for some alignment programs.
    
    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

    Three output files are produced per fasta in the output directory 
    called Output_Translation_Fasta_Files:
    
    	[NAME]_All.fasta - Contains all sequences, pass and fail.
    	
    	[NAME]_Passed.fasta - Contains only sequences that passed translation.
    	
    	[NAME]_Failed.fasta - Contains only sequences that failed translation.

    A summary log file called Log_Sequences_Filtered.txt is written which indicates how many
    sequences passed translation and failed translation for each fasta file processed.

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
import os
import subprocess as sp
import shutil
import argparse
import re
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Coding_Translation_Tests: This function processes an unaligned fasta file. If necessary, it will adjust sequences with 'N's
    so they are divisible by 3 to complete the final codon. Sequences are translated in all
    forward frames to check for the presence of stop codons. If stop codons are detected
    in all frames, it will next check to see if there is only one stop codon and it is present in the
    final two codon positions of the sequence. If the --rc flag is included the translation will also
    be performed for the reverse complement, however if your sequences are all correctly oriented this
    is not recommended. The final sequence is written such that the first base represents the 
    first codon position (based on the detected reading frame) and the sequence is divisible 
    by three to ensure complete codons.

    The translation table should be specified with the --table flag. All NCBI translation table options 
    are available, and can be selected using integers or the shortcut terms provided. If the --table 
    flag is omitted, the default will be to use the Standard translation table. 

    Sequences that have more than one stop codon in all frames fail the translation test. Since no
    reading frame is better than any other, the sequence will be written as is (adjusted for length
    to be divisible by three) to the relevant files, noted by an '*' appended at the end of the sequence
    description line.

    The most important aspects of this script are the adjusting of sequences to the first codon position
    and the completion of the final codon. This is required for some alignment programs. 
    
    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

    Three output files are produced per fasta in the output directory 
    called Output_Translation_Fasta_Files:
    
    	[NAME]_All.fasta - Contains all sequences, pass and fail.
    	
    	[NAME]_Passed.fasta - Contains only sequences that passed translation.
    	
    	[NAME]_Failed.fasta - Contains only sequences that failed translation.

    A summary log file called Log_Sequences_Filtered.txt is written which indicates how many
    sequences passed translation and failed translation for each fasta file processed.
    
    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files. Follow labeling format: NAME.fasta")
    parser.add_argument("--table", default="standard", choices=["standard","vertmtdna","invertmtdna","yeastmtdna","plastid","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"], help="REQUIRED for 'translate' filter: Specifies translation table.")
    parser.add_argument("--rc", action='store_true', help="OPTIONAL: In addition to forward frames, use reverse complement for translation tests. (*Not recommended if sequences have been adjusted to same orientation)")

    return parser.parse_args()

def translation_val(tcode):
    '''
    Convert shortcut translation table terms into the true NCBI
    nomenclature using a dictionary structure. Term needed for 
    translation with biopython.
    '''
    tdict = {"standard":"Standard", "vertmtdna":"Vertebrate Mitochondrial", "invertmtdna":"Invertebrate Mitochondrial", "yeastmtdna":"Yeast Mitochondrial", "plastid":"11", "1":1, "2":2, "3":3, "4":4, "5":5, "6":6, "7":7, "8":8, "9":9, "10":10, "11":11, "12":12, "13":13, "14":14, "15":15, "16":16, "17":17, "18":18, "19":19, "20":20, "21":21, "22":22, "23":23, "24":24, "25":25, "26":26, "27":27, "28":28, "29":29, "30":30, "31":31}
    tval = tdict[tcode]
    return tval

def directory_orf_adjust(in_dir, tcode, rc):
    '''
    Iterates over files in a directory (in_dir) to locate those with
    extension '.fasta'or '.fa' and executes the orf_adjust function
    for each file based on translation table (tcode) and reverse
    complement option (rc). Produces a summary file called
    Sequence_Filtered.txt that contains the fasta name
    and the number of sequences that passed and failed
    translation. Moves all output files to the output
    directory: /Output_Translation_Fasta_Files.
    '''
    print "\n\n--------------------------------------------------------------------------------------"
    print "\nBeginning translation tests and length adjustments.\n"
    print "--------------------------------------------------------------------------------------"
    os.chdir(in_dir)

    out_dir = 'Output_Translation_Fasta_Files'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    drop_name = "Log_Sequences_Filtered.txt"
    with open(drop_name, 'a') as fh_drop_log:
        fh_drop_log.write("Locus\tSeqs_Passed\tSeqs_Failed\n")
    tb = datetime.now()
    f_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])
    for f in f_list:
        summary = orf_adjust(f, tcode, rc)
        with open(drop_name, 'a') as fh_drop_log:
            fh_drop_log.write("{}\t{}\t{}\n".format(summary[0],summary[1],summary[2]))
        output = [o for o in os.listdir('.') if o.endswith('All.fasta') or o.endswith('Passed.fasta') or o.endswith('Failed.fasta')]
        for o in output:
            shutil.move(o, out_dir)
    shutil.move(drop_name, out_dir)
    tf = datetime.now()
    te = tf - tb
    print "\n\n--------------------------------------------------------------------------------------"
    print "\nFinished translation tests and length adjustments."
    print "Total time to process all fasta files: {0} (H:M:S)\n".format(te)
    print "--------------------------------------------------------------------------------------\n\n"

def orf_adjust(fasta_file, tcode, rc):
    '''
    This function processes an unaligned fasta file (fasta_file). If necessary, it will adjust sequences with 'N's
    so they are divisible by 3 to complete the final codon position. Sequences are translated using 
    the specified translation table (tcode) in all forward frames to check for the presence of stop codons. 
    If stop codons are detected in all frames, it will next check to see if there is only one stop codon 
    and it is present in the final two codon positions of the sequence. If the --rc flag is included the 
    translation will also be performed for the reverse complement, however if your sequences are all 
    correctly oriented this is not recommended. The final sequence written such that in the final 
    sequence the first base represents the first codon position (based on the detected reading frame) 
    and the sequence is divisible by three to ensure complete codons.

    Sequences that have more than one stop codon in all frames fail the translation test. Since no
    reading frame is better than any other, the sequence will be written as is (adjusted for length
    to be divisible by three) to the relevant files, noted by an '*' appended at the end of the sequence
    description line.

    Three output files are produced:
    [input fasta name]_All.fasta - Contains all sequences, pass and fail.
    [input fasta name]_Passed.fasta - Contains only sequences that passed translation.
    [input fasta name]_Failed.fasta - Contains only sequences that failed translation.
    '''
    print "\n\nAdjusting reading frames of sequences in {}".format(fasta_file)

    #initiate counts of sequences passing or failing translations
    tpass = int(0)
    tfail = int(0)

    #labels for output files
    names = fasta_file.split('_')
    out_fasta1 = "{}_All.fasta".format(names[0])
    out_fasta2 = "{}_Passed.fasta".format(names[0])
    out_fasta3 = "{}_Failed.fasta".format(names[0])

    #index fasta file and iterate
    fasta_dict = SeqIO.index(fasta_file, "fasta")
    for record in fasta_dict:

        #adjust sequence length to be divisible by three for translation check
        #take all sequences in frame 1 (if possible)
        newseq = fasta_dict[record].seq
        n = Seq("N")
        if len(newseq) % 3 == 0:
            pass
            #print newseq
        elif len(newseq) %3 != 0:
            newseq = newseq + n
            if len(newseq) % 3 == 0:
                pass
                #print newseq
            elif len(newseq) % 3 != 0:
                newseq = newseq + n
                #print newseq
        
        #perform translation with table selected by user
        #translate to all forward frames, constantly adjusting to keep multiples of 3
        frame1 = newseq[:].translate(table=tcode)
        f2 = newseq+n
        frame2 = f2[1:].translate(table=tcode)
        f3 = f2+n
        frame3 = f3[2:].translate(table=tcode)
        #reverse complement and translate to all frames, constantly adjusting to keep multiples of 3
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
        #for each frame - [stop codon count, sequence object in that frame, frame string label, seq slice divisible by 3 based on frame]
        #if rc option is included, use the rc as well
        if rc is False:
            frames_list = [ [f1stops,frame1,"Frame 1",newseq[:]],
                      [f2stops,frame2,"Frame 2",f2[1:]],
                      [f3stops,frame3,"Frame 3",f3[2:]] ]
                
        elif rc is True:
            frames_list = [ [f1stops,frame1,"Frame 1",newseq[:]],
                      [f2stops,frame2,"Frame 2",f2[1:]],
                      [f3stops,frame3,"Frame 3",f3[2:]],
                      [rcf1stops,rcframe1,"RC Frame 1",rcseq[:]],
                      [rcf2stops,rcframe2,"RC Frame 2",rf2[1:]],
                      [rcf3stops,rcframe3,"RC Frame 3",rf3[2:]] ]

        #multi-conditional search statement
        #first, we search all frames to see if they have 0 stop codons, if so we use that sequence frame
        #second, we search all frames to see if a single stop codon occurs in the final two positions, if so we take that frame
        #third, if above two searches fail we take the first (and longest) sequence available
        match = False
        for sublist in frames_list:
            if match is False:        
                #check if frame has no stop codons
                if sublist[0] == int(0):
                    tpass+=1
                    print "\t{0}: Translation {1}".format(fasta_dict[record].id, sublist[2])
                    #make sure that any final codon composed entirely of Ns is removed from seq
                    write_seq = re.sub('NNN$', '', str(sublist[3]))
                    with open(out_fasta1, 'a') as fh_all:
                        fh_all.write(">{0}\n{1}\n".format(fasta_dict[record].description, write_seq) )
                    with open(out_fasta2, 'a') as fh_pass:
                        fh_pass.write(">{0}\n{1}\n".format(fasta_dict[record].description, write_seq) )
                    match = True
        if match is False:
            for sublist in frames_list:
                if match is False:
                    #check if frame has one stop codon and it is in the last two positions of the sequence
                    if sublist[0] == int(1) and '*' in sublist[1][-3:]:
                        tpass+=1
                        print "\t{0}: Translation {1} with single tailing stop codon".format(fasta_dict[record].id, sublist[2])
                        for i in frames_list:
                            print "\t\t\t{} stop codons: {}".format(i[2],i[0])
                        #write_seq = re.sub('NNN$', '', str(sublist[3]))
                        #will remove last two codons (including the detected stop codon)
                        write_seq = str(sublist[3])[:-6]
                        with open(out_fasta1, 'a') as fh_all:
                            fh_all.write( ">{0}\n{1}\n".format(fasta_dict[record].description, write_seq))
                        with open(out_fasta2, 'a') as fh_pass:
                            fh_pass.write( ">{0}\n{1}\n".format(fasta_dict[record].description, write_seq))
                        match = True

        #action if all else fails
        if match is False:
            print "\t{}: WARNING - Failed translation test".format(fasta_dict[record].id)
            for i in frames_list:
                print "\t\t\t{0} stop codons: {1}".format(i[2],i[0])
            print '\t\t\tAdding * to end of the sequence description for {0} in {1}'.format(fasta_dict[record].id,out_fasta1)
            tfail+=1
            write_seq = re.sub('NNN$', '', str(newseq))
            with open(out_fasta1, 'a') as fh_all:
                fh_all.write(">{0}\n{1}\n".format("{}*".format(fasta_dict[record].description), write_seq) )
            with open(out_fasta3, 'a') as fh_fail:
                fh_fail.write(">{0}\n{1}\n".format("{}*".format(fasta_dict[record].description), write_seq) )

    summary = [names[0],tpass,tfail]
    return summary
                        
#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    if args.table is None:
        tcode = "Standard"
    else:
        tcode = translation_val(args.table)

    directory_orf_adjust(args.in_dir, tcode, args.rc)

if __name__ == '__main__':
    main()
