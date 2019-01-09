''''
SuperCRUNCH: Reference_Blast_Extract module

Usage: python Reference_Blast_Extract.py -i [directory containing both files] REQUIRED
										-d [name of reference fasta to make blast database] REQUIRED
										-e [name of empirical fasta to blast to db and trim] REQUIRED
										-b ["blastn", "blastn-short", "dc-megablast", "megablast"] REQUIRED
										--max_hits [integer of maximum number of matches to record] (OPTIONAL)
                                        -m ["span", "nospan", "all"] (OPTIONAL; DEFAULT = span)
										
    Reference_Blast_Extract: Create a blast database and blast an empirical file to the database. The max number of hits
    to the reference allowed is optionally user defined, with the default set to no limit. 
    Query sequences without matches are removed. Query sequences with matches are trimmed 
    to match the resulting blast coordinates. If a single query sequence has multiple 
    coordinates, all overlapping bp coordinates are merged, which may result in a single 
    interval or several non-overlapping intervals. A secondary merging strategy is
    employed because multiple intervals (non-overlapping coordinates) can result from
    two main reasons, including a large stretch or stretches of N's in sequence, or paralogy.
    If paralogy, we want to exclude paralogous intervals, but if N's, we could try to bridge
    those intervals if the number of N's is reasonably low and will allow for good sequence
    alignments downstream. There are three strategies for dealing with multiple intervals:
    
    A) Default option (-m span). Merge intervals if they are <100 bp apart. If this still results
    in multiple non-overlapping intervals (separated by too many N's or from paralogous hits), then
    take the longest interval as the target sequence.
        
    B) 'Single' option (-m nospan). For multiple non-overlapping intervals (separated by too many
    N's or from paralogous hits), then simply take the longest interval as the sequence. This ensures
    that all final sequences are free of N's.
        
    C) 'All' option (-m all). Keeps all intervals and writes the final sequence from these.
    If the non-overlapping intervals result just from stretches of N's, this will
    essentially strip the sequence of N's. However, if it results from paralogous hits,
    then the final sequence will be a composite sequence of target and paralogous sequence.
    Since there is no easy way to objectively check, this is risky and can result in
    spurious sequences. It will be visible in the log file, which lists all coordinates
    used to write a give sequence, and could be useful for finding duplicates or
    paralogy, which often results in sequences much longer than the references.

    The query sequence is then extracted based on the resulting blast interval(s) retained
    after the secondary merge strategy, and written to a new output fasta file.

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
    
    [NAME]_extracted.fasta - The output fasta file containing the extracted sequences.
    
    Log_File_[NAME].txt - A tab delimited text file containing the following columns:
    
    			Accn	Original_Length	Retained_Length	Start_Coordinate	End_Coordinate
                
                The file is populated with data for every record with blast results, such that
                the starting and final sequence lengths are provided along with the interval(s)
                used to extract the final sequence.


    Note - The blast output format 6 is tab-delimited with the following 'hidden' headers:
    qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore


-------------------------
For Python 2.7
Python modules required:
	-blastn (installed in path)
	-BioPython (using SeqIO module)
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
    Reference_Blast_Extract: Create a blast database and blast an empirical file to the database. The max number of hits
    to the reference allowed is optionally user defined, with the default set to no limit. 
    Query sequences without matches are removed. Query sequences with matches are trimmed 
    to match the resulting blast coordinates. If a single query sequence has multiple 
    coordinates, all overlapping bp coordinates are merged, which may result in a single 
    interval or several non-overlapping intervals. A secondary merging strategy is
    employed because multiple intervals (non-overlapping coordinates) can result from
    two main reasons, including a large stretch or stretches of N's in sequence, or paralogy.
    If paralogy, we want to exclude paralogous intervals, but if N's, we could try to bridge
    those intervals if the number of N's is reasonably low and will allow for good sequence
    alignments downstream. There are three strategies for dealing with multiple intervals:
    
    A) Default option (-m span). Merge intervals if they are <100 bp apart. If this still results
    in multiple non-overlapping intervals (separated by too many N's or from paralogous hits), then
    take the longest interval as the target sequence.
        
    B) 'Single' option (-m nospan). For multiple non-overlapping intervals (separated by too many
    N's or from paralogous hits), then simply take the longest interval as the sequence. This ensures
    that all final sequences are free of N's.
        
    C) 'All' option (-m all). Keeps all intervals and writes the final sequence from these.
    If the non-overlapping intervals result just from stretches of N's, this will
    essentially strip the sequence of N's. However, if it results from paralogous hits,
    then the final sequence will be a composite sequence of target and paralogous sequence.
    Since there is no easy way to objectively check, this is risky and can result in
    spurious sequences. It will be visible in the log file, which lists all coordinates
    used to write a give sequence, and could be useful for finding duplicates or
    paralogy, which often results in sequences much longer than the references.

    The query sequence is then extracted based on the resulting blast interval(s) retained
    after the secondary merge strategy, and written to a new output fasta file.

    The type of blast search is defined by the user, the most popular being "dc-megablast",
    "blastn",and "megablast". Note that blastn searches include word size of 11, whereas 
    megablast searches include a word size of 28, making blastn more appropriate for inter-
    species searches and megablast more appropriate for closely related or intraspecific 
    searches. The discontiguous megablast is better at producing non-fragmented hits for 
    divergent sequences using similar word sizes as blastn, and is preferable to blastn.

    Empirical fasta file should be labeled as 'NAME.fasta', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

	Outputs Files: [NAME]_blast_output.txt - Contains blast search results in output format 6;
    [NAME]_extracted.fasta - The output fasta file containing the filtered sequences.
    Log_File_[NAME].txt - A tab delimited text file containing the following columns:
    Accn	Original_Length	Retained_Length	Start_Coordinate	End_Coordinate
    The file is populated with data for every record with blast results, such that
    the starting and final sequence lengths are provided along with the interval(s)
    used to extract the final sequence.
    
    DEPENDENCIES: Python: BioPython; Executables in path: blast tools (makeblastdb, blastn).
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory containing the reference and the empirical fasta files.")
    parser.add_argument("-d", "--blast_db", required=True, help="REQUIRED: The name of the reference fasta file that will be used to create the blast database. File name only (rather than full path), should be located in the input directory.")
    parser.add_argument("-e", "--emp_fasta", required=True, help="REQUIRED: The name of an empirical fasta to blast to the database to prune sequences. File name only (rather than full path), should be located in the input directory. Follow labeling format: NAME.fasta")
    parser.add_argument("-b", "--blast_task", required=True, choices=["blastn", "blastn-short", "dc-megablast", "megablast"], help="REQUIRED: The type of blast search to conduct.")
    parser.add_argument("--max_hits", type=int, default=None, help="OPTIONAL: The maximum number of blast matches allowed per input sequence. May want to set < 300.")
    parser.add_argument("-m", "--merge_strategy", required=False, default="span", choices=["span", "nospan", "all"], help="DEFAULT = 'span'. The strategy for dealing with multiple non-overlapping blast coordinates. See documentation for details.")

    return parser.parse_args()

#-----------------------------------------------------------------------------------------------

def make_db(blast_db):
    print "\n---------------------------------------------------------------------------\n"
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    print '\n\n', mdb_str
    proc = sp.call(mdb_str, shell=True)

def blastn_to_db(task, blast_db, emp_fasta, outname, max_hits):
    print "\n---------------------------------------------------------------------------\n"
    if max_hits is None:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 > {3}".format(task, blast_db, emp_fasta, outname)
    	print '\n\n', blast_str
    else:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} > {4}".format(task, blast_db, emp_fasta, max_hits, outname)
    	print '\n\n', blast_str
    proc = sp.call(blast_str, shell=True)

def merge_coords(interval_list):
    '''
    Takes a list of sublists where each sublist is composed of an
    interval [start integer, stop integer]. Assumes that for
    intervals the first number is lower than the second number.
    The sublists are sorted by the start integer, low to
    high, and the search begins with the first interval/sublist.
    If the intervals overlap the starting interval list is
    updated with the new maximum stop integer, and the search
    continues. If there is no overlap, the non-overlapping
    interval becomes the new target for the search, which
    ends after the sublists have been iterated over.
    
    Example:
    This list of sublists (sorted):
    [ [1,6], [5,15], [7,12], [20,25], [22,35], [35,40] ]
    
    Creates the merged list:
    [ [1, 6, 15, 15], [20, 25, 35, 40] ]

    Which is cleaned to create the final list:
    [ [1, 15], [20, 40] ]
    '''
    #sort the list of lists by their first element (start coordinate)
    interval_list.sort(key=lambda x: x[0])
    #initiate empty list to store updated coordinates
    merged = []
    #iterate over lists
    for interval in interval_list:
        #if merged list is empty, add first interval/sublist
        #otherwise, see if the stop coord of the next interval/sublist
        #is less than the start coord of this (merged) interval
        if not merged or merged[-1][-1] < interval[0]:
            #if there is no overlap, add the non-overlapping interval
            #to the merged list to begin the next round of searching
            merged.append(interval)
        #if there is overlap
        else:
            #find the highest stop coord between the current (merged)
            #interval and the overlapping interval
            int_end = max(merged[-1][-1], interval[-1])
            #append that number to the end of the merged list entry,
            #growing the (merged) interval by the highest stop integer
            #found so far
            merged[-1].append(int_end)
    #take first and last element from each sublist in merged
    #this produces one start and stop coord for every interval
    cleaned = []
    for m in merged:
        cleaned.append([m[0], m[-1]])

    return cleaned

def span_coords(coords):
    '''
    Check if multiple bp intervals are within 100bp of each
    other, and if so, merge them and return a new list
    of updated coordinates.
    '''
    new_coords = []
    for i in range(0, (len(coords) - 1)):
        #check if the end coordinate of the first entry is within 100bp
        #of the start coordinate of the subsequent entry
        #if so, merge
        if (int(coords[i+1][0]) - int(coords[i][1])) <= int(100):
            updated = [coords[i][0],coords[i+1][1]]
            new_coords.append(updated)
        else:
            new_coords.append(coords[i])
            new_coords.append(coords[i+1])
            
    new_coords.sort(key=lambda x: x[0])
    
    return new_coords

def choose_coord(coords):
    '''
    If a coordinate list has multiple non-overlapping intervals,
    find the longest interval and return the single coordinate
    set as a list.
    '''
    diff_index = []
    for i, j in enumerate(coords):
        diff = j[-1]-j[0]
        diff_index.append([i,diff])
        
    diff_index.sort(key=lambda x: x[1], reverse=True)
    selected = coords[(diff_index[0][0])]
    return selected

def one_coord_span(merged_coords):
    '''
    Function for finding the most sensible blast hit to use.
    If multiple non-overlapping coordinates are present,
    will first attempt to connect intervals that are separated
    by less than 100bp. If multiple non-overlapping coordinates
    persist, then a long stretch (or stretches) of N's is present in the
    sequence (making it low-quality) or it resulted from
    paralogous hits. In this case, the longest interval is
    retained (which should represent the target locus).
    '''
    #NOTE: merged_coords always returns a list of sublists,
    #even if only one sublist, ie [[start,stop]]
    #must return same data structure here
    if len(merged_coords) >= int(2):
        spanned = span_coords(merged_coords)
        if len(spanned) >= int(2):
            spanned = merge_coords(spanned)
                           
        if len(spanned) >= int(2):
            final_coord = choose_coord(spanned)
        else:
            final_coord = spanned[0]
    else:
        final_coord = merged_coords[0]

    return final_coord

def one_coord_nospan(merged_coords):
    '''
    Function for finding the most sensible blast hit to use.
    If multiple non-overlapping coordinates are present,
    then a long stretch (or stretches) of N's is present in the
    sequence (making it low-quality) or it resulted from
    paralogous hits. In this case, the longest interval is
    retained (which should represent the target locus). Skips the
    attempt to span sequences (could result in N's in final seq).
    '''
    if len(merged_coords) >= int(2):
        final_coord = choose_coord(merged_coords)
    else:
        final_coord = merged_coords[0]

    return final_coord

def get_accs_and_contents(f):
    '''
    Extract information from blast output file.
    '''
    #First create a set of accn numbers, which we will turn into
    #a list, sort, and iterate over. At the same time, store the file
    #contents as a split and stripped list of items.
    print "Parsing contents in {}...\n".format(f)
    accn_set = set()
    file_contents = []
    with open(f, 'r') as fh_in:
        for line in fh_in:
            #list comprehension to strip whitespace from all items results from splitting 
            #by the tab character, resulting in a list of file line contents
            parts = [l.strip() for l in line.split('\t')]
            accn_set.add(parts[0])
            file_contents.append(parts)
    accn_list = list(accn_set)
    accn_list.sort()

    return file_contents, accn_list

def parse_blastn_output6_NoSelfHits(outname, merge_strategy):
    '''
    Get information from blast output file.
    '''
    print "\n---------------------------------------------------------------------------\n"
    #extract file lines and sorted list of accession numbers from blast output file
    file_contents, accn_list = get_accs_and_contents(outname)
    
    #create a list containing only the important info: [accn, [start1 bp, stop1 bp]]
    #if multiple sections blasted to reference then list will look like this: [accn, [start1 bp, stop1 bp], [start2 bp, stop2 bp] ...]
    parsing_list = []
    #initiate counter for accns:
    count = int(0)
    #look up for every accn in the list
    for a in accn_list:
        #test if number of accns processed divisible by 100, if so print number
        #ie an on-screen progress report
        if count != int(0) and count % int(100) == 0:
        	print "\t\t\tProcessed hits for {} accessions...".format(count)
        #initiate list to store all bp coordinates found for this accn 
        coord_list = []
        #iterate over the file contents
        for l in file_contents:
            #find accn match, but exclude self-hits if they exist
            if l[0] == a and l[1] != a:
                #add coords to list based on the index of the split line
                #subtract 1 to adjust indexing for output
                coord1 = int(l[6]) - int(1)
                coord2 = int(l[7]) - int(1)
                coord_list.append([coord1, coord2])
        #use merge function to merge all overlapping coordinates found in blast output file
        merged_list = merge_coords(coord_list)
        
        #Initiate list to store information for accession number:
        #Under -m options span or nospan - [accn number, [coordinate set1]]
        #Under -m option all - [accn number, [coordinate set1], [coordinateset2], ...]
        summary_list = [a]

        #This can return single intervals like: [1,500]
        #But can also return multiple intervals like: [[0, 445], [509, 1088]]

        #Multiple intervals (non-overlapping coordinates) can result from two main reasons,
        #including a large stretch or stretches of N's in sequence, or paralogy.
        #If paralogy, we want to exclude paralogous intervals, but if N's, we could try
        #to bridge those intervals if the number of N's is reasonably low and will allow
        #for good sequence alignments downstream.
        
        #I've come up with three strategies for dealing with multiple intervals:
        #1) Default option (-m span) Merge intervals if they are <100 bp apart. If this still results
        #in multiple non-overlapping intervals (separated by too many N's or from paralogous hits), then
        #take the longest interval as the target sequence.
        
        #2) 'Pure' option (-m nospan) For multiple non-overlapping intervals (separated by too many
        #N's or from paralogous hits), then simply take the longest interval as the sequence. This ensures
        #that all final sequences are free of N's.
        
        #3) 'All' option (-m all) Keeps all intervals and writes the final sequence from these.
        #If the non-overlapping intervals result just from stretches of N's, this will
        #essentially strip the sequence of N's. However, if it results from paralogous hits,
        #then the final sequence will be a composite sequence of target and paralogous sequence.
        #Since there is no easy way to objectively check, this is risky and can result in
        #spurious sequences. It will be visible in the log file, which lists all coordinates
        #used to write a give sequence, and could be useful for finding duplicates or
        #paralogy, which often results in sequences much longer than the references.
        
        if merge_strategy == "span":
            final_coord = one_coord_span(merged_list)
            summary_list.append(final_coord)

        elif merge_strategy == "nospan":
            final_coord = one_coord_nospan(merged_list)
            summary_list.append(final_coord)

        elif merge_strategy == "all":
            #if accn had more than one non-overlapping merged coord, we add all.
            for m in merged_list:
                summary_list.append(m)
        parsing_list.append(summary_list)
        #add to accn count
        count += 1
    print "Finished retrieving blast coordinates for {} accessions.\n".format(count)
    
    return parsing_list

def pull_records(emp_fasta, parsing_list):
    print "\n---------------------------------------------------------------------------\n"
    print "Writing sequences with updated coordinates...\n"
    #load the fasta file using the indexing function
    #will create a dictionary with accn as keys
    record_dict = SeqIO.index(emp_fasta, "fasta")
    #initiate output files
    autoname = "{}_extracted.fasta".format(emp_fasta.split('.')[0])
    logname = "Log_File_{}.txt".format(emp_fasta.split('.')[0])
    with open(logname, 'a') as fh_log:
        fh_log.write("Accn\tOriginal_Length\tRetained_Length\tCoordinates_Used\n")
    with open(autoname, 'a') as fh_out:
        #initiate count for records sliced and written
        count = int(0)
        #iterate over list with accns and coordinates
        for item in parsing_list:
            #test if number of accns processed divisible by 100, if so print number
            #ie an on-screen progress report
            if count != int(0) and count % int(100) == 0:
                print "\t\t\tFinished writing {} updated records...".format(count)
            #look up the sequence using the accn number as the dictionary key
            fullseq = record_dict[item[0]].seq
            seq_parts = []
            #iterate over parsing list where coordinate lists would start (index 1)
            #if only one list, fine, but if multiple we need to examine them to
            #pull out the appropriate parts of the sequence of interest
            for coord in item[1:]:
                #trim sequence if necessary, convert seq object to string
                #which allows for concatenation later (if multiple intervals present)
                seqslice = str(fullseq[coord[0]:(coord[1] + int(1))])
                seq_parts.append(seqslice)
            #if there are multiple sequence sections, join them here
            #if only one item in list will not throw an error
            newseq = "".join(seq_parts)
            #write information to log file
            with open(logname, 'a') as fh_log:
                fh_log.write("{0}\t{1}\t{2}\t{3}\n".format(item[0], len(fullseq), len(newseq), item[1:]))
            #write to updated fasta file
            fh_out.write(">{}\n{}\n".format(record_dict[item[0]].description, newseq))
            #add to counter
            count += 1
    print "Wrote a total of {0} sequences to {1}\n\n".format(count, autoname) 

#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    os.chdir(args.in_dir)
    tb = datetime.now()
    make_db(args.blast_db)
    outname = "{}_blast_output.txt".format((args.emp_fasta).split('.')[0])
    blastn_to_db(args.blast_task, args.blast_db, args.emp_fasta, outname, args.max_hits)
    parsing_list = parse_blastn_output6_NoSelfHits(outname,args.merge_strategy)
    pull_records(args.emp_fasta, parsing_list)
    tf = datetime.now()
    te = tf - tb
    print "Total time to process {0}: {1} (H:M:S)\n\n".format(args.emp_fasta,te)

if __name__ == '__main__':
    main()
