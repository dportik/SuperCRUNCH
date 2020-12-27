''''
SuperCRUNCH: Contamination_Filter module

    Contamination_Filter: Create a blast database from 'contamination' sequences 
    and blast an empirical file to this database. Any query sequence that scores 
    a pident of > 95% with length > 100 bp is considered to be a match to a 
    'contamination' sequence. These bad sequences are written to a separate fasta file, 
    whereas all other passing sequences are written to a new filtered fasta file.

    For example, use human mtDNA as the contamination seq database and blast 
    sequences of a different group of organisms (reptiles, amphibians, birds, 
    etc) to the database to look for human contaminated sequences. 

    The type of blast search is defined by the user. For the purpose of finding 
    contaminated sequences, megablast may be preferable.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
   
    Note - The blast output format 6 is tab-delimited with the following 
           'hidden' headers:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore


-------------------------
Compatible with Python 2.7 & 3.7
Dependencies:
	-blast+ (installed in path)
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
import shutil
import os
import subprocess as sp
import argparse
import operator
from datetime import datetime
from Bio import SeqIO


def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Contamination_Filter: Create a blast database from 'contamination' sequences 
    and blast an empirical file to this database. Any query sequence that scores 
    a pident of > 95% with length > 100 bp is considered to be a match to a 
    'contamination' sequence. These bad sequences are written to a separate fasta file, 
    whereas all other passing sequences are written to a new filtered fasta file.

    For example, use human mtDNA as the contamination seq database and blast 
    sequences of a different group of organisms (reptiles, amphibians, birds, 
    etc) to the database to look for human contaminated sequences. 

    The type of blast search is defined by the user. For the purpose of finding 
    contaminated sequences, megablast may be preferable.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
    
    DEPENDENCIES: Python: BioPython; Executables in path: blast+ tools 
    (makeblastdb, blastn).

	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory containing the "
                            "reference and the empirical fasta files")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to write "
                            "output files.")
    
    parser.add_argument("-e", "--emp_fasta",
                            required=False,
                            help="OPTION 1: The name of an empirical fasta to blast to "
                            "the database to identify 'contaminated' sequences. File name "
                            "only (rather than full path), should be located in the input "
                            "directory. Follow labeling format: NAME.fasta")
    
    parser.add_argument("-d", "--blast_db",
                            required=False,
                            help="OPTION 1: The name of the reference fasta file (of "
                            "sequences from the contamination source) that "
                            "will be used to create the blast database. File name only "
                            "(rather than full path), should be located in the input directory.")
    
    parser.add_argument("--multisearch",
                            required=False,
                            help="OPTION 2: If multiple empirical fasta + blast database "
                            "pairs are present in the -i (input) directory, provide the full path "
                            "to a tab-delimited text file containing the pair information "
                            "to allow all of them to run sequentially. See documentation.")
    
    parser.add_argument("-b", "--blast_task",
                            required=True,
                            choices=["blastn", "blastn-short", "dc-megablast", "megablast"],
                            help="REQUIRED: The type of blast search to conduct.")
        
    parser.add_argument("--max_hits",
                            type=int,
                            default=None,
                            help="OPTIONAL: The maximum number of blast matches allowed "
                            "per input sequence.")
    
    return parser.parse_args()

def make_db(blast_db):
    """
    Generates call string for makeblastd and then
    runs using sp.call() with shell option. 
    """
    print("\n---------------------------------------------------------------------------\n")
    
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    
    print("\n\n{}".format(mdb_str))
    proc = sp.call(mdb_str, shell=True)

def blastn_to_db(task, blast_db, emp_fasta, outname, max_hits):
    """
    Generates call string for blastn for file (f) using the 
    supplied arguments. Runs using sp.call() with shell 
    option. Returns the name of the output file generated by 
    the program.
    """
    print("\n---------------------------------------------------------------------------\n")
    
    if max_hits is None:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 > {3}"
                         .format(task, blast_db, emp_fasta, outname))
    else:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} > {4}"
                         .format(task, blast_db, emp_fasta, max_hits, outname))
            
    print("\n\n{}".format(blast_str))
    proc = sp.call(blast_str, shell=True)


def parse_blastn_output6(outname):
    """
    Get information from blast output file.
    """
    print("\n---------------------------------------------------------------------------\n")
    
    #Store the file contents of the blast output file
    #as a split and stripped list.
    print("Parsing contents in {}...\n".format(outname))
    file_contents = []
    
    with open(outname, 'r') as fh_in:
        for line in fh_in:
            temp_list = []
            #list comprehension to strip whitespace from all items results from splitting 
            #by the tab character, resulting in a list of file line contents
            parts = [l.strip() for l in line.split('\t')]
            file_contents.append(parts)
            
    #sort blast output sublists by reverse order of 'pident'
    file_contents.sort(key=operator.itemgetter(2), reverse=True)

    #create accession lists for contaminated and good seqs
    bad_seqs = []
    #iterate over sorted file contents list
    for l in file_contents:
        #test if pident is > 95% and length is > 100bp (matching contaminated seq)
        if float(l[2]) > 95.0 and int(l[3]) > 100:
            bad_seqs.append(l[0])
            
    return bad_seqs

def pull_records_from_accn_lists(emp_fasta, bad_seqs):
    """
    Function to write sequences to new fasta file based 
    on whether they passed the contam filter or not.
    """
    #load the fasta file using the indexing function
    #will create a dictionary with accn as keys
    record_dict = SeqIO.index(emp_fasta, "fasta")
    
    prefix = emp_fasta.split('.')[0]
    
    #initiate counts for records written
    goodcount, badcount = int(0), int(0)

    print("\n---------------------------------------------------------------------------")
    print("\nWriting sequences that passed contamination filter...\n")

    #first write sequences that passed filter
    goodfasta = "{}_contfiltered.fasta".format(prefix)
    with open(goodfasta, 'a') as fh:
        #iterate over seq records
        for acc in record_dict:
            #check if accession number is not present in bad seq list
            if acc not in bad_seqs:
                #test if number of accns processed is divisible by 100, if so print number (ie progress report for big files)
                if goodcount % int(100) == 0:
                    print("\t\t\tFinished writing {:,} records...".format(goodcount))
                #write record to fasta
                fh.write(">{}\n{}\n".format(record_dict[acc].description, record_dict[acc].seq))
                #add to counter
                goodcount += 1
            
    print("\n---------------------------------------------------------------------------")
    print("\nWriting contaminated sequences...\n")

    #second write any sequences that failed the filter
    badfasta = "{}_contaminated.fasta".format(prefix)
    with open(badfasta, 'a') as fh:
        if len(bad_seqs) >= int(1):
            #iterate over list with accns and coordinates
            for acc in bad_seqs:
                #write record to fasta
                fh.write(">{}\n{}\n".format(record_dict[acc].description, record_dict[acc].seq))
                #add to counter
                badcount += 1
            print("\t\t\tFinished writing {} records...".format(badcount))
            
    #log for bad record IDs
    logname = "{}_Contaminated_Seq_Summary.txt".format(prefix)
    with open(logname, 'a') as fh_out:
        fh_out.write("Wrote a total of {0} sequences to {1}\n\n".format(goodcount, goodfasta))
        fh_out.write("Wrote a total of {0} sequences to {1}, including:\n\n".format(badcount, badfasta))
        if len(bad_seqs) >= int(1):
            for acc in bad_seqs:
                #write record to fasta
                fh_out.write("{}\n".format(record_dict[acc].description))

    print("\n\nWrote a total of {0:,} sequences to {1}".format(goodcount, goodfasta))
    print("Wrote a total of {0:,} sequences to {1}".format(badcount, badfasta))

def cleanup(blastdir, cdir, fdir):
    """
    Moves relevant output files to their output directories. 
    """
    print("\nCleaning up output files...")
    
    [os.remove(f) for f in os.listdir('.') if f.endswith('contaminated.fasta') and os.stat(f).st_size == 0]
    [shutil.move(f, cdir) for f in os.listdir('.') if f.endswith("contaminated.fasta")]
    [shutil.move(f, blastdir) for f in os.listdir('.') if f.endswith(("blast_results.txt", "Contaminated_Seq_Summary.txt"))]
    [shutil.move(f, fdir) for f in os.listdir('.') if f.endswith("contfiltered.fasta")]
    [os.remove(f) for f in os.listdir('.') if f.endswith((".nsq", ".nhr", ".nin", ".ndb", ".not", ".ntf", ".nto"))]
    
    print("\tDone!\n")

def make_dirs(outdir):
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    os.chdir(outdir)
    curpath = os.getcwd()
        
    blastdir = os.path.join(curpath, "Blast-Output-Files")
    if not os.path.exists(blastdir):
        os.mkdir(blastdir)
        
    cdir = os.path.join(curpath, "Contaminated-Sequences")
    if not os.path.exists(cdir):
        os.mkdir(cdir)
        
    fdir = os.path.join(curpath, "Filtered-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)

    return blastdir, cdir, fdir

def main():
    args = get_args()
    tb = datetime.now()

    blastdir, cdir, fdir = make_dirs(args.outdir)
    
    os.chdir(args.indir)
    '''
    make_db(args.blast_db)
    
    outname = "{}_blast_results.txt".format((args.emp_fasta).split('.')[0])
    blastn_to_db(args.blast_task, args.blast_db, args.emp_fasta, outname, args.max_hits)
    
    bad_seqs = parse_blastn_output6(outname)
    pull_records_from_accn_lists(args.emp_fasta, bad_seqs)
    
    cleanup(blastdir, cdir, fdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal time to process {0}: {1} (H:M:S)\n".format(args.emp_fasta, te))
    print("--------------------------------------------------------------------------------------\n\n")
    '''
    if not args.multisearch:
        print("\n\n======================================================================================")
        print("\nProcessing fasta file: {}\n".format(args.emp_fasta.split('.')[0]))
        print("======================================================================================\n\n")
        make_db(args.blast_db)
        outname = "{}_blast_results.txt".format((args.emp_fasta).split('.')[0])
        blastn_to_db(args.blast_task, args.blast_db, args.emp_fasta, outname, args.max_hits)
    
        bad_seqs = parse_blastn_output6(outname)
        pull_records_from_accn_lists(args.emp_fasta, bad_seqs)

        cleanup(blastdir, cdir, fdir)
        tf = datetime.now()
        te = tf - tb
        print("\n\n--------------------------------------------------------------------------------------")
        print("\nTotal time to process {0}: {1} (H:M:S)\n".format(args.emp_fasta, te))
        print("--------------------------------------------------------------------------------------\n\n")    

        
    elif args.multisearch:
        with open(args.multisearch, 'r') as fh:
            pairs = [l.strip().split('\t') for l in fh if l.strip()]
        if pairs:
            for p in pairs:
                print("\n\n======================================================================================")
                print("\nProcessing fasta file: {}\n".format(p[0]))
                print("======================================================================================\n\n")
                b = datetime.now()
                make_db(p[1])
                outname = "{}_blast_results.txt".format((p[0]).split('.')[0])
                blastn_to_db(args.blast_task, p[1], p[0], outname, args.max_hits)

                bad_seqs = parse_blastn_output6(outname)
                pull_records_from_accn_lists(p[0], bad_seqs)

                cleanup(blastdir, cdir, fdir)
                f = datetime.now()
                e = f - b
                print("\nElapsed time: {0} (H:M:S)\n\n".format(e))
                
        else:
            print("\n\nSomething went wrong while attempting to parse file: {}".format(args.multisearch.split('/')[-1]))
            print("Please ensure this is a tab-delimited file with empirical fasta name in column 1 and "
                      "the database fasta name in column 2, and all files are within the -i directory.\n\n")
    
        tf = datetime.now()
        te = tf - tb
        print("\n\n--------------------------------------------------------------------------------------")
        print("\nTotal time to process {0} file(s): {1} (H:M:S)\n".format(len(pairs), te))
        print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()
