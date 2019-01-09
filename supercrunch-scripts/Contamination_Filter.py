''''
SuperCRUNCH: Contamination_Filter module

Usage: python Contamination_Filter.py -i [directory containing both files] REQUIRED
                                    -d [name of reference fasta to make blast database] REQUIRED
                                    -e [name of empirical fasta to blast to db] REQUIRED
                                    -b ["blastn", "blastn-short", "dc-megablast", "megablast"] REQUIRED
                                    --max_hits [integer of maximum number of matches to record] (OPTIONAL)

    Contamination_Filter: Create a blast database from 'contamination' sequences and blast an empirical file 
    to this database. Any query sequence that scores a pident of > 95% with length > 100 bp
    is considered to be a match to the contaminated seqs. These bad sequences are written to the 
    [NAME]_contaminated.fasta file, whereas all other sequences are written to the
    [NAME]_filtered.fasta file. The record IDs of all contaminated sequences are
    written to Log_Contaminated_Seqs.txt.

    For example, use human mtDNA as the contamination seq database and blast sequences of a different
    group of organisms (reptiles, amphibians, birds, etc) to the database to look for human
    contaminated sequences. 

    The type of blast search is defined by the user, the most popular being "dc-megablast",
    "blastn",and "megablast". Note that blastn searches include word size of 11, whereas 
    megablast searches include a word size of 28, making blastn more appropriate for inter-
    species searches and megablast more appropriate for closely related or intraspecific 
    searches. The discontiguous megablast is better at producing non-fragmented hits for 
    divergent sequences using similar word sizes as blastn, and is preferable to blastn.
    For the purpose of finding contaminated sequences, megablast may be preferable.

    Empirical fasta file should be labeled as 'NAME.fasta', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.
   
    Outputs Files:
    
    [NAME]_blast_output.txt - Contains blast search results in output format 6
    
    [NAME]_filtered.fasta - The output fasta file containing the sequences passing the contamination filter.
    
    [NAME]_contaminated.fasta - The output fasta file containing the putatively contaminated sequences.
    
    Log_Contaminated_Seqs.txt - A tab delimited text file containing the record IDs for all bad seqs.
    
    Note - The blast output format 6 is tab-delimited with the following 'hidden' headers:
    qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore


-------------------------
For Python 2.7
Dependencies:
	-blastn (installed in path)
	-BioPython (using SeqIO module here)
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
January 2019
Distributed under the 
GNU General Public Lincense
'''
import sys
import os
import subprocess as sp
import argparse
import operator
from datetime import datetime
from Bio import SeqIO


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Contamination_Filter: Create a blast database from 'contamination' sequences and blast an empirical file 
    to this database. Any query sequence that scores a pident of > 95% with length > 100 bp
    is considered to be a match to the contaminated seqs. These bad sequences are written to the 
    [NAME]_contaminated.fasta file, whereas all other sequences are written to the
    [NAME]_filtered.fasta file. The record IDs of all contaminated sequences are
    written to Log_Contaminated_Seqs.txt.

    For example, use human mtDNA as the contamination seq database and blast sequences of a different
    group of organisms (reptiles, amphibians, birds, etc) to the database to look for human
    contaminated sequences. 

    The type of blast search is defined by the user, the most popular being "dc-megablast",
    "blastn",and "megablast". Note that blastn searches include word size of 11, whereas 
    megablast searches include a word size of 28, making blastn more appropriate for inter-
    species searches and megablast more appropriate for closely related or intraspecific 
    searches. The discontiguous megablast is better at producing non-fragmented hits for 
    divergent sequences using similar word sizes as blastn, and is preferable to blastn.
    
    Empirical fasta file should be labeled as 'NAME.fasta', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.
   
    Outputs Files:
    
    [NAME]_blast_output.txt - Contains blast search results in output format 6
    
    [NAME]_filtered.fasta - The output fasta file containing the sequences passing the contamination filter.
    
    [NAME]_contaminated.fasta - The output fasta file containing the putatively contaminated sequences.
    
    Log_Contaminated_Seqs.txt - A tab delimited text file containing the record IDs for all bad seqs.
    
    DEPENDENCIES: Python: BioPython; Executables in path: blast tools (makeblastdb, blastn).

	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory containing the reference and the empirical fasta files")
    parser.add_argument("-d", "--blast_db", required=True, help="REQUIRED: The name of the reference fasta file that will be used to create the blast database. File name only (rather than full path), should be located in the input directory.")
    parser.add_argument("-e", "--emp_fasta", required=True, help="REQUIRED: The name of an empirical fasta to blast to the database to prune sequences. File name only (rather than full path), should be located in the input directory. Follow labeling format: NAME.fasta")
    parser.add_argument("-b", "--blast_task", required=True, choices=["blastn", "blastn-short", "dc-megablast", "megablast"], help="REQUIRED: The type of blast search to conduct")
    parser.add_argument("--max_hits", type=int, default=None, help="OPTIONAL: The maximum number of blast matches allowed per input sequence.")
    return parser.parse_args()

#-----------------------------------------------------------------------------------------------

def make_db(blast_db):
    '''
    Function to make a blast database from file 'blast_db'.
    '''
    print "\n---------------------------------------------------------------------------\n"
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    print '\n\n', mdb_str
    proc = sp.call(mdb_str, shell=True)

def blastn_to_db(task, blast_db, emp_fasta, outname, max_hits):
    '''
    Function blast a fasta file (emp_fast) to a database (blast_db) using a
    particular blast algorithm (task) with or without a maximum number of hits (max_hits).
    '''    
    print "\n---------------------------------------------------------------------------\n"
    if max_hits is None:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 > {3}".format(task, blast_db, emp_fasta, outname)
    	print '\n\n', blast_str
    else:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} > {4}".format(task, blast_db, emp_fasta, max_hits, outname)
    	print '\n\n', blast_str
    proc = sp.call(blast_str, shell=True)


def parse_blastn_output6(outname):
    '''
    Get information from blast output file (outname).
    '''
    print "\n---------------------------------------------------------------------------\n"
    
    #Store the file contents of the blast output file
    #as a split and stripped list.
    print "Parsing contents in {}...\n".format(outname)
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
    #load the fasta file using the indexing function
    #will create a dictionary with accn as keys
    record_dict = SeqIO.index(emp_fasta, "fasta")

    #initiate counts for records written
    goodcount = int(0)
    badcount = int(0)

    #write good seqs
    print "\n---------------------------------------------------------------------------\nWriting sequences that passed contamination filter...\n"
    goodfasta = "{}_filtered.fasta".format(emp_fasta.split('.')[0])
    with open(goodfasta, 'a') as fh_out:
        #iterate over seq records
        for acc in record_dict:
            #check if accession number is not present in bad seq list
            if acc not in bad_seqs:
                #test if number of accns processed is divisible by 100, if so print number (ie progress report for big files)
                if goodcount % int(100) == 0:
                    print "\t\t\tFinished writing {} records...".format(goodcount)
                #write record to fasta
                fh_out.write(">{}\n{}\n".format(record_dict[acc].description, record_dict[acc].seq))
                #add to counter
                goodcount += 1
            
    #write bad seqs, if they occur, to fasta file
    print "\n---------------------------------------------------------------------------\nWriting contaminated sequences...\n"
    badfasta = "{}_contaminated.fasta".format(emp_fasta.split('.')[0])
    with open(badfasta, 'a') as fh_out:
        if len(bad_seqs) >= int(1):
            #iterate over list with accns and coordinates
            for acc in bad_seqs:
                #write record to fasta
                fh_out.write(">{}\n{}\n".format(record_dict[acc].description, record_dict[acc].seq))
                #add to counter
                badcount += 1
            print "\t\t\tFinished writing {} records...".format(badcount)
            
    #log for bad record IDs
    logname = "Log_Contaminated_Seqs.txt"
    with open(logname, 'a') as fh_out:
        fh_out.write("Wrote a total of {0} sequences to {1}\n\n".format(goodcount, goodfasta))
        fh_out.write("Wrote a total of {0} sequences to {1}, including:\n\n".format(badcount, badfasta))
        if len(bad_seqs) >= int(1):
            for acc in bad_seqs:
                #write record to fasta
                fh_out.write("{}\n".format(record_dict[acc].description))

    print "\n\nWrote a total of {0} sequences to {1}".format(goodcount, goodfasta) 
    print "Wrote a total of {0} sequences to {1}\n\n".format(badcount, badfasta) 

#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    os.chdir(args.in_dir)
    tb = datetime.now()
    make_db(args.blast_db)
    outname = "{}_blast_output.txt".format((args.emp_fasta).split('.')[0])
    blastn_to_db(args.blast_task, args.blast_db, args.emp_fasta, outname, args.max_hits)
    bad_seqs = parse_blastn_output6(outname)
    pull_records_from_accn_lists(args.emp_fasta, bad_seqs)
    tf = datetime.now()
    te = tf - tb
    print "Total time to process {0}: {1} (H:M:S)\n\n".format(args.emp_fasta,te)

if __name__ == '__main__':
    main()
