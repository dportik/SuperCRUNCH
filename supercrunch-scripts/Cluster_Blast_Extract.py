''''
SuperCRUNCH: Cluster_Blast_Extract module

    Cluster_Blast_Extract: A stringent similarity filter for filtering sequences.
    Each fasta file should contain putative sequences for a SINGLE locus. For
    each fasta file, the sequences are first clustered using cd-hit-est. A blast
    database is constructed from the largest sequence cluster (which is assumed
    to contain the true sequences for that locus), and all sequences are blasted
    to the sequence database (including those in the cluster used to generate the database).
    The type of blast search is defined by the user, the most popular being 
    "dc-megablast", "blastn",and "megablast". Note that blastn searches include a 
    word size of 11, whereas megablast searches include a word size of 28, making 
    blastn more appropriate for inter-species searches and megablast more appropriate 
    for closely related or intraspecific searches. The discontiguous megablast is 
    better at producing non-fragmented hits for divergent sequences using similar 
    word sizes as blastn, and is preferable to blastn. The number of threads to use
    for blast searches can be specified using the --threads flag. For all hits of a 
    particular sequence (exlcuding self-hits), the blast coordinates are collected 
    and merged. If multiple non-overlapping intervals are found (ex. 10-40, 45-80), 
    the merging strategy selected by the user (-m) is used (see documentation for 
    details, including explanation of --bp_brige argument). Finally, the query 
    sequence is extracted based on the blast coordinates and written to a new 
    filtered fasta file. Sequences that do not produce blast hits fail the 
    similarity filter and are written to a separate fasta file. The output files
    from each step are written to separate directories which will appear in the
    main output directory specified.

    Main Output Files:
    
    [NAME]_blast_output.txt - Contains blast search results in output format 6
    
    [NAME]_extracted.fasta - The output fasta file containing the extracted sequences.
    
    Log_File_[NAME].txt - A tab delimited text file containing the following columns:
    
    			Accn	Original_Length	Retained_Length	Start_Coordinate	End_Coordinate
                
                The file is populated with data for every record with blast results, such that
                the starting and final sequence lengths are provided along with the interval(s)
                used to extract the final sequence.

    DEPENDENCIES: Python: BioPython; Executables in path: blast+ tools 
    (makeblastdb, blastn), cd-hit-est. ***Note the regular version of cd-hit-est
    produced the same results across multiple runs, but the openmp version
    compiled from source gave inconsistent results. It should probably be avoided
    until it is fixed (e.g., compile using 'make openmp=no'). 

-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
Dependencies:
    -blast+ (installed in path)
    -cd-hit-est (installed in path)
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
import itertools
import operator
from datetime import datetime
from Bio import SeqIO

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Cluster_Blast_Extract: A stringent similarity filter for filtering sequences.
    Each fasta file should contain putative sequences for a SINGLE locus. For
    each fasta file, the sequences are first clustered using cd-hit-est. A blast
    database is constructed from the largest sequence cluster (which is assumed
    to contain the true sequences for that locus), and all sequences are blasted
    to the sequence database (including those in the cluster used to generate the database).
    The type of blast search is defined by the user, the most popular being 
    "dc-megablast", "blastn",and "megablast". Note that blastn searches include a 
    word size of 11, whereas megablast searches include a word size of 28, making 
    blastn more appropriate for inter-species searches and megablast more appropriate 
    for closely related or intraspecific searches. The discontiguous megablast is 
    better at producing non-fragmented hits for divergent sequences using similar 
    word sizes as blastn, and is preferable to blastn. The number of threads to use
    for blast searches can be specified using the --threads flag. For all hits of a 
    particular sequence (excluding self-hits), the blast coordinates are collected 
    and merged. If multiple non-overlapping intervals are found (ex. 10-40, 50-80), 
    the merging strategy selected by the user (-m) is used (see documentation for 
    details, including explanation of --bp_brige argument). Finally, the query 
    sequence is extracted based on the blast coordinates and written to a new 
    filtered fasta file. Sequences that do not produce blast hits fail the 
    similarity filter and are written to a separate fasta file. The output files
    from each step are written to separate directories which will appear in the
    main output directory specified.

    DEPENDENCIES: Python: BioPython; Executables in path: blast+ toolkit 
    (makeblastdb, blastn), cd-hit-est. ***Note the regular version of cd-hit-est
    produced the same results across multiple runs, but the openmp version
    compiled from source gave inconsistent results. It should probably be avoided
    until it is fixed (e.g., compile using 'make openmp=no'). 
	-----------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--in_dir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the parsed, locus-specific fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="OPTIONAL: The full path to an existing directory to write "
                            "output files.")
    
    parser.add_argument("-b", "--blast_task",
                            required=True,
                            choices=["blastn", "blastn-short", "dc-megablast", "megablast"],
                            help="REQUIRED: The type of blast search to conduct.")
    
    parser.add_argument("--max_hits",
                            type=int,
                            default=None,
                            help="OPTIONAL: The maximum number of blast matches allowed "
                            "per input sequence. (May want to set < 300)")
    
    parser.add_argument("-m", "--merge_strategy",
                            required=False,
                            default="span", choices=["span", "nospan", "all"],
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

def cluster(f):
    """
    Generates call string for cd-hit-est and then
    runs using sp.call() with shell option. Returns
    name of the output file generated by the program.
    """
    name = f.split('.')[0]
    
    call_string = "cd-hit-est -i {0}.fasta -o {0}_clustering_out -c 0.8 -n 4 -M 16000 -d 0".format(name)
    proc = sp.call(call_string, shell=True)

    #produces file called [name]_clustering_out.clstr
    cluster_file = "{0}_clustering_out.clstr".format(name)
    
    return cluster_file

def delimit(iterable, splitstring):
    """
    Function to split cd-hit-est cluster files in a sensible way.
    Groups the items in the list using a keyword, and creates a 
    list of sublists that are grouped accordingly. Returns the
    list of sublists.
    """
    return [list(g) for k,g in itertools.groupby(iterable,lambda x:x in splitstring) if not k]

def parse_clusters(f, cluster_file, fdict):
    """
    Reads a cluster file and generates a list of the contents, where
    the string Cluster denotes the beginning of a new sequence cluster.
    The delimit() function is used to create a list of sublists based on the 
    appearance of the string Cluster in the original list. The list of sublists
    is reverse sorted by length to identify the largest cluster. A fasta file
    for this cluster is then generated by finding the accessions of the contained
    sequences and using the fasta dictionary (fdict) to write the sequences
    with biopython. The name of the cluster fasta file is returned.    
    """
    print("\n--------------------------------------------------------------------------------------\n")
    print('Parsing clusters for {}.'.format(cluster_file))
    lines = []

    #A clstr file delimits clusters with a >, then numbers all included accession numbers
    #>Cluster 0
    #0	616nt, >JN881132.1... at +/86.20%
    #1	561nt, >KU765220.1... at +/96.08%
    #2	558nt, >KU765300.1... at +/96.59%
    #3	522nt, >JF818216.1... at *
    #Collect the file contents then delimit lists based on presence of 'Cluster' in text
    with open(cluster_file, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                lines.append("Cluster")
            else:
                lines.append(line.strip())
    #use delimit function to parse clusters into a list of lists
    sublists = delimit(lines,("Cluster",))

    #retrieve the largest cluster which will be used as
    #blast databse
    sublists.sort(key=len, reverse=True)
    print('\tFound {} cluster(s).'.format(len(sublists)))
    print('\tLargest cluster contains {:,} sequences.'.format(len(sublists[0])))

    #write largest cluster to new fasta file
    cluster_db = "{0}_clusterDB.fasta".format(f.split('.')[0])
    
    with open(cluster_db, 'a') as fh:
        #access the biggest cluster, get all seqs
        for string in sublists[0]:
            if '>' in string and '...' in string:
                accession = string.split('>')[1].split('...')[0]
                #use fdict to look up accession and write record
                fh.write((fdict[accession]).format("fasta"))
                
    return cluster_db
            
def make_bdb(blast_db):
    """
    Generates call string for makeblastd and then
    runs using sp.call() with shell option. 
    """
    print("\n--------------------------------------------------------------------------------------\n")
    
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    
    print(mdb_str)
    proc = sp.call(mdb_str, shell=True)
    
def blastn_to_db(f, blast_db, task, max_hits, threads):
    """
    Generates call string for blastn for file (f) using the 
    supplied arguments. Runs using sp.call() with shell 
    option. Returns the name of the output file generated by 
    the program.
    """
    outname = "{}_blast_results.txt".format(f.split('.')[0])
    
    if threads is None:
        threads = 1
        
    if max_hits is None:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 -num_threads {3} > {4}"
                         .format(task, blast_db, f, threads, outname))
    else:
    	blast_str = ("blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} -num_threads {4} > {5}"
                         .format(task, blast_db, f, max_hits, threads, outname))
        
    print("\n\n{}\n".format(blast_str))
    proc = sp.call(blast_str, shell=True)
    
    return outname

def get_accs_and_contents(f):
    """
    Extract information from blast output file (f),
    including a list of split lines and a list of 
    accession numbers. Returns both lists.
    """
    print("Parsing contents in {}...".format(f))
    accn_set = set()
    file_contents = []
    
    with open(f, 'r') as fh_in:
        for line in fh_in:
            parts = [l.strip() for l in line.split('\t')]
            accn_set.add(parts[0])
            file_contents.append(parts)
            
    accn_list = sorted(accn_set)
    
    return file_contents, accn_list

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
    interval_list.sort(key=lambda x: x[0])
    merged = []
    
    for interval in interval_list:
        if not merged or merged[-1][-1] < interval[0]:
            merged.append(interval)
            
        else:
            int_end = max(merged[-1][-1], interval[-1])
            merged[-1].append(int_end)
            
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
        #if so, merge. X can be set by user (arg --bp_bridge) or default is 100.
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
        diff_index.append([i, diff])
        
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

def parse_blastn_output6_NoSelfHits(blast_output, merge_strategy, bp_bridge):
    """
    Get information from blast output file based on merge_strategy
    selected. See annotations below.
    """
    print("\n--------------------------------------------------------------------------------------\n")
    #extract file lines and sorted list of accession numbers from blast output file
    file_contents, accn_list = get_accs_and_contents(blast_output)
    
    #create a list containing only the important info: [accn, [start1 bp, stop1 bp]]
    #if multiple sections blasted to reference then list will look like this:
    #[accn, [start1 bp, stop1 bp], [start2 bp, stop2 bp] ...]
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
            #find accn match, but exclude self-hits
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
            
    print("\tFinished retrieving blast coordinates for {:,} accessions.".format(count))
    
    return parsing_list

def pull_records(f, fdict, parsing_list):
    """
    Function to write sequences to new fasta file based 
    on the blast coordinates retrieved for those sequences.
    """
    print("\nWriting sequences with updated coordinates...")
    
    #create set of accession numbers from empirical fasta
    emp_accs = set(list(fdict.keys()))
    #initiate empty set to populate with accession for seqs extracted
    extracted_accs = set()
    #initiate output files
    autoname = "{}_extracted.fasta".format(f.split('.')[0])
    logname = "Log_File_{}.txt".format(f.split('.')[0])
    badname = "Log_BadSeqs_{}.fasta".format(f.split('.')[0])
    
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
            fullseq = fdict[item[0]].seq
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
            fh_out.write(">{}\n{}\n".format(fdict[item[0]].description, newseq))
            #add to counter
            count += 1
    #write file of sequences that failed
    badseqs = emp_accs - extracted_accs
    if len(badseqs) >= 1:
        with open(badname, 'a') as fh_out:
            for acc in badseqs:
                fh_out.write(">{}\n{}\n".format(fdict[acc].description, (fdict[acc].seq)))

    print("\tWrote a total of {0:,} sequences to {1}.".format(count, autoname))
    
    if len(badseqs) >= 1:
        print("\n{0:,} starting sequence(s) did not pass orthology filtering and are written to {1}."
                  .format(len(badseqs), badname))
        
    elif len(badseqs) == int(0):
        print("\nAll starting sequences passed filtering!")


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
    clustdir = os.path.join(curpath, "Clustering-Output-Files")
    if not os.path.exists(clustdir):
        os.mkdir(clustdir)
        
    blastdir = os.path.join(curpath, "Blast-Output-Files")
    if not os.path.exists(blastdir):
        os.mkdir(blastdir)
        
    logdir = os.path.join(curpath, "Blast-Trimming-Log-Files")
    if not os.path.exists(logdir):
        os.mkdir(logdir)
        
    fdir = os.path.join(curpath, "Filtered-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)

    return clustdir, blastdir, logdir, fdir

def cleanup(cluster_outs, clustdir, blastdir, logdir, fdir):
    """
    Moves relevant output files to their output directories. 
    """
    print("\n\nCleaning up output files...")

    [shutil.move(f, clustdir) for f in cluster_outs]
    [shutil.move(f, logdir) for f in os.listdir('.') if f.startswith(("Log_File", "Log_BadSeqs"))]
    [shutil.move(f, blastdir) for f in os.listdir('.') if f.endswith("blast_results.txt")]
    [os.remove(f) for f in os.listdir('.') if os.stat(f).st_size == 0 and f.endswith("_extracted.fasta")]
    [shutil.move(f, fdir) for f in os.listdir('.') if f.endswith("_extracted.fasta")]
    [os.remove(f) for f in os.listdir('.') if f.endswith((".nsq", ".nhr", ".nin"))]
    
    print("\tDone!\n")


def cbe_runner(f, blast_task, max_hits, merge_strategy,
                   threads, bp_bridge, clustdir, blastdir, logdir, fdir):
    """
    Function to run the major steps for a particular file (f).
    This includes:
    1) load the sequences to memory in dictionary
       format using biopython
    2) cluster sequences using cd-hit-est with
       the cluster() function
    3) parse the clusters to identify the largest cluster
       using the parse_clusters() function
    4) create a blast database from the largest cluster 
       using the make_bdb() function
    5) blast all sequences to the database using the
       blastn_to_db() function
    6) obtain blast coordinates for all sequences with hits
       using the parse_blastn_output6_NoSelfHits() function
    7) write the filtered sequences to a new fasta (using their coordinates)
       with the pull_records() function
    8) clean up temporary files, moving all outputs to relevant
       output directories using the cleanup() function
    """
    b = datetime.now()
    print("\n\n======================================================================================")
    print("\nProcessing fasta file: {}\n".format(f))
    print("======================================================================================\n\n")
    
    fdict = SeqIO.index(f, "fasta")
    
    cluster_file = cluster(f)
    
    cluster_db = parse_clusters(f, cluster_file, fdict)
    
    make_bdb(cluster_db)
    
    blast_output = blastn_to_db(f, cluster_db, blast_task, max_hits, threads)
    
    parsing_list = parse_blastn_output6_NoSelfHits(blast_output, merge_strategy, bp_bridge)
    
    if parsing_list:
        pull_records(f, fdict, parsing_list)

    cluster_outs = [cluster_file, cluster_db, cluster_file.replace(".clstr", "")]
    
    cleanup(cluster_outs, clustdir, blastdir, logdir, fdir)
    
    f = datetime.now()
    e = f - b
    print("\nElapsed time: {0} (H:M:S)\n\n\n".format(e))
    

def main():
    args = get_args()
    tb = datetime.now()

    clustdir, blastdir, logdir, fdir = make_dirs(args.outdir)
    
    os.chdir(args.in_dir)
    fastas = sorted([f for f in os.listdir('.') if f.endswith((".fa", ".fasta"))])

    for f in fastas:
        cbe_runner(f, args.blast_task, args.max_hits, args.merge_strategy,
                       args.threads, args.bp_bridge, clustdir, blastdir, logdir, fdir)
    
    clean_fastas = sorted([f for f in os.listdir(fdir) if f.endswith('.fasta')])
    os.chdir(args.outdir)
    with open("Filtered_Fasta_File_List.txt", 'a') as fh:
        for f in clean_fastas:
            fh.write("{}\n".format(f))

    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")
    
if __name__ == '__main__':
    main()
