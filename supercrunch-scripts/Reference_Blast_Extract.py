''''
SuperCRUNCH: Reference_Blast_Extract module
										
    Reference_Blast_Extract: A stringent similarity filter for filtering sequences.
    Creates a blast database from user-selected reference sequences and blasts 
    a supplied fasta file to this database. 

    The type of blast search is defined by the user, the most popular being 
    "dc-megablast", "blastn",and "megablast". Note that blastn searches include a 
    word size of 11, whereas megablast searches include a word size of 28, making 
    blastn more appropriate for inter-species searches and megablast more appropriate 
    for closely related or intraspecific searches. The discontiguous megablast is 
    better at producing non-fragmented hits for divergent sequences using similar 
    word sizes as blastn, and is preferable to blastn. The number of threads to use
    for blast searches can be specified using the --threads flag. For all hits of a 
    particular sequence (excluding self-hits), the blast coordinates are collected 
    and merged. If multiple non-overlapping intervals are found (ex. 10-40, 45-80), 
    the merging strategy selected by the user (-m) is used (see documentation for 
    details, including explanation of --bp_brige argument). Finally, the query 
    sequence is extracted based on the blast coordinates and written to a new 
    filtered fasta file. Sequences that do not produce blast hits fail the 
    similarity filter and are written to a separate fasta file. The output files
    from each step are written to separate directories which will appear in the
    main output directory specified.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
   
    Note - The blast output format 6 is tab-delimited with the following 'hidden' headers:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
Dependencies:
    -blast+ (installed in path)
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
import shutil
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
    Reference_Blast_Extract: A stringent similarity filter for filtering sequences.
    Creates a blast database from user-selected reference sequences and blasts 
    a supplied fasta file to this database. 

    The type of blast search is defined by the user, the most popular being 
    "dc-megablast", "blastn",and "megablast". Note that blastn searches include a 
    word size of 11, whereas megablast searches include a word size of 28, making 
    blastn more appropriate for inter-species searches and megablast more appropriate 
    for closely related or intraspecific searches. The discontiguous megablast is 
    better at producing non-fragmented hits for divergent sequences using similar 
    word sizes as blastn, and is preferable to blastn. The number of threads to use
    for blast searches can be specified using the --threads flag. For all hits of a 
    particular sequence (excluding self-hits), the blast coordinates are collected 
    and merged. If multiple non-overlapping intervals are found (ex. 10-40, 45-80), 
    the merging strategy selected by the user (-m) is used (see documentation for 
    details, including explanation of --bp_brige argument). Finally, the query 
    sequence is extracted based on the blast coordinates and written to a new 
    filtered fasta file. Sequences that do not produce blast hits fail the 
    similarity filter and are written to a separate fasta file. The output files
    from each step are written to separate directories which will appear in the
    main output directory specified.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
    
    DEPENDENCIES: Python: BioPython; Executables in path: blast+ tools (makeblastdb, blastn).
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory containing "
                            "the reference(s) and the empirical fasta file(s).")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to write "
                            "output files.")
    
    parser.add_argument("-e", "--emp_fasta",
                            required=False,
                            help="OPTION 1: The name of an empirical fasta to blast to "
                            "the database to prune sequences. File name only (rather than "
                            "full path), should be located in the input directory. Follow "
                            "labeling format: NAME.fasta")
    
    parser.add_argument("-d", "--blast_db",
                            required=False,
                            help="OPTION 1: The name of the reference fasta file that "
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
                            "per input sequence. May want to set < 300.")
    
    parser.add_argument("-m", "--merge_strategy",
                            required=False,
                            default="span",
                            choices=["span", "nospan", "all"],
                            help="DEFAULT = 'span'. The strategy for dealing with multiple "
                            "non-overlapping blast coordinates. See documentation for details.")
    
    parser.add_argument("--bp_bridge",
                            required=False,
                            default="100",
                            type=int,
                            help="OPTIONAL: The base pair distance used to merge non-overlapping "
                            "blast coordinates for -m span (default = 100).")
    
    parser.add_argument("--threads",
                            default=None,
                            help="OPTIONAL: Specifies number of threads to use for blastn "
                            "(recommended: 4).")

    return parser.parse_args()

#-----------------------------------------------------------------------------------------------

def make_db(blast_db):
    """
    Generates call string for makeblastd and then
    runs using sp.call() with shell option. 
    """
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    
    print("\n\n{}".format(mdb_str))
    proc = sp.call(mdb_str, shell=True)

def blastn_to_db(task, blast_db, emp_fasta, outname, max_hits, threads):
    """
    Generates call string for blastn for file (f) using the 
    supplied arguments. Runs using sp.call() with shell 
    option. Returns the name of the output file generated by 
    the program.
    """
    print("\n---------------------------------------------------------------------------\n")

    if threads is None:
        threads = 1
    
    if max_hits is None:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 -num_threads {3} > {4}"
                         .format(task, blast_db, emp_fasta, threads, outname))
    else:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} -num_threads {4} > {5}"
                         .format(task, blast_db, emp_fasta, max_hits, threads, outname))
            
    print("\n\n{}".format(blast_str))
    proc = sp.call(blast_str, shell=True)

def merge_coords(interval_list):
    """
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
    """
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

def span_coords(coords, bp_bridge):
    """
    Check if multiple bp intervals are within X bp of each
    other, and if so, merge them and return a new list
    of updated coordinates. X is determined by argument
    bp_bridge, which defaults to 100 unless supplied by user.
    """
    new_coords = []
    for i in range(0, (len(coords) - 1)):
        #check if the end coordinate of the first entry is within X bp
        #of the start coordinate of the subsequent entry
        #if so, merge. X can be set by user (arg --bp_bridge; default is 100)
        if (int(coords[i+1][0]) - int(coords[i][1])) <= int(bp_bridge):
            updated = [coords[i][0],coords[i+1][1]]
            new_coords.append(updated)
            
        else:
            new_coords.append(coords[i])
            new_coords.append(coords[i+1])
            
    new_coords.sort(key=lambda x: x[0])
    
    return new_coords

def choose_coord(coords):
    """
    If a coordinate list has multiple non-overlapping intervals,
    find the longest interval and return the single coordinate
    set as a list.
    """
    diff_index = []
    for i, j in enumerate(coords):
        diff = j[-1]-j[0]
        diff_index.append([i,diff])
        
    diff_index.sort(key=lambda x: x[1], reverse=True)
    selected = coords[(diff_index[0][0])]
    
    return selected

def one_coord_span(merged_coords, bp_bridge):
    """
    Function for finding the most sensible blast coordinates to use.
    If multiple non-overlapping coordinates are present,
    will first attempt to connect intervals that are separated
    by less than X bp (default 100). If multiple intervals
    persist, then a long stretch of N's is present in the
    sequence (making it low-quality) or it resulted from
    distinct hits from one sequence. In this case, the longest 
    interval is retained (which should represent the target locus).
    """
    #NOTE: merged_coords always returns a list of sublists,
    #even if only one sublist, ie [[start,stop]]
    #must return same data structure here
    if len(merged_coords) >= int(2):
        spanned = span_coords(merged_coords, bp_bridge)
        
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
    """
    Function for finding the most sensible blast coordinates to use.
    If multiple non-overlapping coordinates are present,
    then a long stretch (or stretches) of N's is present in the
    sequence (making it low-quality) or it resulted from
    paralogous hits. In this case, the longest interval is
    retained (which should represent the target locus). Skips the
    attempt to span sequences (could result in N's in final seq).
    """
    if len(merged_coords) >= int(2):
        final_coord = choose_coord(merged_coords)
        
    else:
        final_coord = merged_coords[0]

    return final_coord

def get_accs_and_contents(f):
    """
    Extract information from blast output file (f),
    including a list of split lines and a list of 
    accession numbers. Returns both lists.
    """
    #First create a set of accn numbers, which we will turn into
    #a list, sort, and iterate over. At the same time, store the file
    #contents as a split and stripped list of items.
    print("Parsing contents in {}...\n".format(f))
    accn_set = set()
    file_contents = []
    
    with open(f, 'r') as fh_in:
        for line in fh_in:
            #list comprehension to strip whitespace from all items results from splitting 
            #by the tab character, resulting in a list of file line contents
            parts = [l.strip() for l in line.split('\t')]
            accn_set.add(parts[0])
            file_contents.append(parts)
            
    accn_list = sorted(accn_set)

    return file_contents, accn_list

def parse_blastn_output6_NoSelfHits(outname, merge_strategy, bp_bridge):
    """
    Get information from blast output file based on merge_strategy
    selected. See annotations below.
    """
    print("\n---------------------------------------------------------------------------\n")
    #extract file lines and sorted list of accession numbers from blast output file
    file_contents, accn_list = get_accs_and_contents(outname)
    
    #create a list containing only the important info: [accn, [start1 bp, stop1 bp]]
    #if multiple sections blasted to reference then list will look like this:
    # [accn, [start1 bp, stop1 bp], [start2 bp, stop2 bp] ...]
    parsing_list = []
    #initiate counter for accns:
    count = int(0)
    #look up for every accn in the list
    for a in accn_list:
        #test if number of accns processed divisible by 100, if so print number
        #ie an on-screen progress report
        if count != int(0) and count % int(100) == 0:
        	print("\t\tProcessed hits for {:,} accessions...".format(count))
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
                
        #make sure coord list not empty (ie blast results were pure self hits)
        if coord_list:
            #use merge function to merge all overlapping coordinates found in blast output file
            merged_list = merge_coords(coord_list)

            #Initiate list to store information for accession number:
            summary_list = [a]

            #deal with coordinate merging strategy
            if merge_strategy == "span":
                final_coord = one_coord_span(merged_list, bp_bridge)
                summary_list.append(final_coord)

            elif merge_strategy == "nospan":
                final_coord = one_coord_nospan(merged_list)
                summary_list.append(final_coord)

            elif merge_strategy == "all":
                #if accn had more than one non-overlapping merged coord, we add all.
                for m in merged_list:
                    summary_list.append(m)
                    
            parsing_list.append(summary_list)
            count += 1
    print("\nFinished retrieving blast coordinates for {:,} accessions.\n".format(count))
    
    return parsing_list

def pull_records(emp_fasta, parsing_list):
    """
    Function to write sequences to new fasta file based 
    on the blast coordinates retrieved for those sequences.
    """
    print("\n---------------------------------------------------------------------------\n")
    print("Writing sequences with updated coordinates...\n")
    #load the fasta file using the indexing function
    #will create a dictionary with accn as keys
    record_dict = SeqIO.index(emp_fasta, "fasta")
    #create set of accession numbers from empirical fasta
    emp_accs = set(list(record_dict.keys()))
    #initiate empty set to populate with accession for seqs extracted
    extracted_accs = set()
    #initiate output files
    autoname = "{}_extracted.fasta".format(emp_fasta.split('.')[0])
    logname = "Log_File_{}.txt".format(emp_fasta.split('.')[0])
    badname = "Log_BadSeqs_{}.fasta".format(emp_fasta.split('.')[0])
    
    with open(logname, 'a') as fh_log:
        fh_log.write("Accn\tOriginal_Length\tRetained_Length\tCoordinates_Used\n")
    with open(autoname, 'a') as fh_out:
        #initiate count for records sliced and written
        count = int(0)
        #iterate over list with accns and coordinates
        for item in parsing_list:
            extracted_accs.add(item[0])
            #test if number of accns processed divisible by 100, if so print number
            #ie an on-screen progress report
            if count != int(0) and count % int(100) == 0:
                print("\t\tFinished writing {:,} updated records...".format(count))
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
                fh_log.write("{0}\t{1}\t{2}\t{3}\n"
                                 .format(item[0], len(fullseq), len(newseq), item[1:]))
            #write to updated fasta file
            fh_out.write(">{}\n{}\n".format(record_dict[item[0]].description, newseq))
            #add to counter
            count += 1
            
    #write file of sequences that failed
    badseqs = emp_accs - extracted_accs
    if len(badseqs) >= 1:
        with open(badname, 'a') as fh_out:
            for acc in badseqs:
                fh_out.write(">{}\n{}\n".format(record_dict[acc].description, (record_dict[acc].seq)))

    print("\nWrote a total of {0:,} sequences to {1}.".format(count, autoname))
    
    if len(badseqs) >= 1:
        print("{0:,} starting sequence(s) did not pass similarity filtering and are written to {1}.\n\n"
                  .format(len(badseqs), badname))
        
    elif len(badseqs) == int(0):
        print("All starting sequences passed similarity filtering!\n\n")
        
def cleanup(blastdir, logdir, fdir):
    """
    Moves relevant output files to their output directories. 
    """
    print("\nCleaning up output files...")

    [shutil.move(f, logdir) for f in os.listdir('.') if f.startswith(("Log_File", "Log_BadSeqs"))]
    [shutil.move(f, blastdir) for f in os.listdir('.') if f.endswith("blast_results.txt")]
    [shutil.move(f, fdir) for f in os.listdir('.') if f.endswith("_extracted.fasta")]
    [os.remove(f) for f in os.listdir('.') if f.endswith((".nsq", ".nhr", ".nin"))]
    
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
        
    logdir = os.path.join(curpath, "Blast-Trimming-Log-Files")
    if not os.path.exists(logdir):
        os.mkdir(logdir)
        
    fdir = os.path.join(curpath, "Filtered-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)

    return blastdir, logdir, fdir

def main():
    tb = datetime.now()
    args = get_args()
    
    blastdir, logdir, fdir = make_dirs(args.outdir)
    
    os.chdir(args.indir)
    
    if not args.multisearch:
        print("\n\n======================================================================================")
        print("\nProcessing fasta file: {}\n".format(args.emp_fasta.split('.')[0]))
        print("======================================================================================\n\n")
        make_db(args.blast_db)
        outname = "{}_blast_results.txt".format((args.emp_fasta).split('.')[0])
        blastn_to_db(args.blast_task, args.blast_db, args.emp_fasta, outname, args.max_hits, args.threads)
    
        parsing_list = parse_blastn_output6_NoSelfHits(outname, args.merge_strategy, args.bp_bridge)
        pull_records(args.emp_fasta, parsing_list)

        cleanup(blastdir, logdir, fdir)
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
                blastn_to_db(args.blast_task, p[1], p[0], outname, args.max_hits, args.threads)

                parsing_list = parse_blastn_output6_NoSelfHits(outname, args.merge_strategy, args.bp_bridge)
                pull_records(p[0], parsing_list)

                cleanup(blastdir, logdir, fdir)
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
