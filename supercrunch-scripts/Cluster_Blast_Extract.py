'''
SuperCRUNCH: Cluster_Blast_Extract module

Usage: python Cluster_Blast_Extract.py  -i [directory containing fasta files] REQUIRED
                                        -b ["blastn", "blastn-short", "dc-megablast", "megablast"] REQUIRED
                                        -m ["span", "nospan", "all"] (OPTIONAL; DEFAULT = span)
                                        --max_hits [integer] (OPTIONAL)
									
    Cluster_Blast_Extract: A stringent paralogy filter for sorting GenBank data of a putative locus.
    
    Fasta files should contain putative sequences for a SINGLE locus. Each fasta file therefore
    represents a different gene/locus of interest.
    To be recognized, a fasta file must be labeled as "NAME.fasta". The NAME portion 
    should not contain any periods or spaces, but can contain underscores. 

										
    Major steps performed for each fasta file found:

    1. Perform cd-hit-est to create sequence clusters based on similarity.

    2. Create a fasta file for each cluster recovered.

    3. Sort fasta clusters by number of records contained (high to low).

    4. Create a blast database from the cluster containing the most records (which presumably
       contains homologous sequences for that locus).
   
    5. Perform blast searches of all clusters to this database. This includes the cluster used 
       to create the blast database, but self-hits are ignored. The max number of hits
       to the reference allowed is optionally user defined, with the default set to no limit.
   
       The type of blast search is defined by the user, the most popular being "dc-megablast",
       "blastn",and "megablast". Note that blastn searches include word size of 11, whereas 
       megablast searches include a word size of 28, making blastn more appropriate for inter-
       species searches and megablast more appropriate for closely related or intraspecific 
       searches. The discontiguous megablast is better at producing non-fragmented hits for 
       divergent sequences using similar word sizes as blastn, and is preferable to blastn.
       This should ideally be used for inter- specific searches, replacing blastn.
   
    6. Merge all blast results, and for every accession/record with blast results merge all
       overlapping coordinates (excluding self hits). If multiple non-overlapping intervals
       are found (ex. 10-40, 45-80), a secondary merging strategy is employed because multiple
       intervals (non-overlapping coordinates) can result from two main reasons, including a
       large stretch or stretches of N's in sequence, or paralogy. If paralogy, we want to
       exclude paralogous intervals, but if N's, we could try to bridge those intervals if
       the number of N's is reasonably low and will allow for good sequence alignments
       downstream. There are three strategies for dealing with multiple intervals:
    
        A) Default option (-m span). Merge intervals if they are <100 bp apart. If this still
        results in multiple non-overlapping intervals (separated by too many N's or from
        paralogous hits), then take the longest interval as the target sequence.
        
        B) 'Single' option (-m nospan). For multiple non-overlapping intervals (separated by
        too many N's or from paralogous hits), then simply take the longest interval as the
        sequence. This ensures that all final sequences are free of N's.
        
        C) 'All' option (-m all). Keeps all intervals and writes the final sequence from these.
        If the non-overlapping intervals result just from stretches of N's, this will
        essentially strip the sequence of N's. However, if it results from paralogous hits,
        then the final sequence will be a composite sequence of target and paralogous sequence.
        Since there is no easy way to objectively check, this is risky and can result in
        spurious sequences. It will be visible in the log file, which lists all coordinates
        used to write a give sequence, and could be useful for finding duplicates or
        paralogy, which often results in sequences much longer than the references.
   
    7. Create a new filtered fasta file for each sequence with results. The query sequence is
       extracted based on the resulting blast interval(s) retained after the secondary merge
       strategy, and written to a new output fasta file. Sequences that did not result in
       significant blast hits did not pass the paralogy filter and are excluded. Checking 
       these I have found they are often paralogous, mislabeled, contaminated (ie wrong 
       organism), or the wrong type (ex. mtDNA or mRNA). Clusters that failed to blast 
       to the database will have a [fasta name]_Out_Cluster_#_blastn_results.txt output 
       file with zero size (empty file).  You can examine these clusters and perform a 
       blast search to a global GenBank database to validate these observations.

    Several output directories are created in the main directory which contain different output files.


    01_Clustering_Results:

    	[fasta name].clstr - The messy output file from cd-hit-est that contains cluster information.



    02_Parsed_Results:

    	[fasta name]_Out_Cluster_# - Fasta files representing each numbered cluster. These are 
    								used for blast searches, and the largest fasta is used to
	    							create the blast database.
								

    03_Blast_Results:

	    [fasta name]_blast_results_merged.txt - The complete, merged blastn results for all the
                								clusters found for that locus.
								
	    [fasta name]_Out_Cluster_#_blast_results.txt - The blastn results for all the numbered
    			                					cluster fasta files found for that locus.


    04_Trimmed_Results:

	    [fasta name]_extracted.fasta - A new fasta file that is created from the sequences extracted
		    						based on the blastn coordinates found. Only sequences with 
			    					significant blast results will be included, and depending
				    				on the blast coordinates found, sequences may or may not 
					    			have been trimmed from their original length (see file below).
								
    	Log_File_[fasta name].txt - A summary file that shows all records included in the filtered
	    							fasta file. Includes accession number, the original length
		    						of the unfiltered sequence, the final length of the extracted
			    					sequence, and the bp coordinates used to extract the sequence
				    				(which come directly from the merged blast hits for that 
					    			record). The headers of this file include:
						    		Accn	Original_Length		Retained_Length		Coordinates_Used

-------------------------
For Python 2.7
Python modules required:
	-BioPython (using SeqIO module)
Other dependencies:
	-cd-hit-est (installed in path)
	-blastn (installed in path)
	    -must be able to call makeblastdb from command line
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
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Cluster_Blast_Extract: A stringent paralogy filter for sorting GenBank data of a putative locus.
    
    Fasta files should contain putative sequences for a SINGLE locus. Each fasta file therefore
    represents a different gene/locus of interest.
    To be recognized, a fasta file must be labeled as "NAME.fasta". The NAME portion 
    should not contain any periods or spaces, but can contain underscores. 

										
    Major steps performed for each fasta file found:

    1. Perform cd-hit-est to create sequence clusters based on similarity.

    2. Create a fasta file for each cluster recovered.

    3. Sort fasta clusters by number of records contained (high to low).

    4. Create a blast database from the cluster containing the most records (which presumably
       contains homologous sequences for that locus).
   
    5. Perform blast searches of all clusters to this database. This includes the cluster used 
       to create the blast database, but self-hits are ignored. The max number of hits
       to the reference allowed is optionally user defined, with the default set to no limit.
   
       The type of blast search is defined by the user, the most popular being "dc-megablast",
       "blastn",and "megablast". Note that blastn searches include word size of 11, whereas 
       megablast searches include a word size of 28, making blastn more appropriate for inter-
       species searches and megablast more appropriate for closely related or intraspecific 
       searches. The discontiguous megablast is better at producing non-fragmented hits for 
       divergent sequences using similar word sizes as blastn, and is preferable to blastn.
       This should ideally be used for inter- specific searches, replacing blastn.
   
    6. Merge all blast results, and for every accession/record with blast results merge all
       overlapping coordinates (excluding self hits). If multiple non-overlapping intervals
       are found (ex. 10-40, 45-80), a secondary merging strategy is employed because multiple
       intervals (non-overlapping coordinates) can result from two main reasons, including a
       large stretch or stretches of N's in sequence, or paralogy. If paralogy, we want to
       exclude paralogous intervals, but if N's, we could try to bridge those intervals if
       the number of N's is reasonably low and will allow for good sequence alignments
       downstream. There are three strategies for dealing with multiple intervals:
    
        A) Default option (-m span). Merge intervals if they are <100 bp apart. If this still
        results in multiple non-overlapping intervals (separated by too many N's or from
        paralogous hits), then take the longest interval as the target sequence.
        
        B) 'Single' option (-m nospan). For multiple non-overlapping intervals (separated by
        too many N's or from paralogous hits), then simply take the longest interval as the
        sequence. This ensures that all final sequences are free of N's.
        
        C) 'All' option (-m all). Keeps all intervals and writes the final sequence from these.
        If the non-overlapping intervals result just from stretches of N's, this will
        essentially strip the sequence of N's. However, if it results from paralogous hits,
        then the final sequence will be a composite sequence of target and paralogous sequence.
        Since there is no easy way to objectively check, this is risky and can result in
        spurious sequences.
   
    7. Create a new filtered fasta file for each sequence with results. The query sequence is
       extracted based on the resulting blast interval(s) retained after the secondary merge
       strategy, and written to a new output fasta file. Sequences that did not result in
       significant blast hits did not pass the paralogy filter and are thus excluded.

    DEPENDENCIES: Python: BioPython; Executables in path: blast tools (makeblastdb, blastn), 
    cd-hit-est (***NOT openmp compiled version). The clustering results were unable to be 
    replicated with the openmp version of cd-hit-est, but were consistent and replicable using
    the regular version (compiled using 'make openmp=no').
	-----------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the parsed, locus-specific fasta files.")
    parser.add_argument("-b", "--blast_task", required=True, choices=["blastn", "blastn-short", "dc-megablast", "megablast"], help="REQUIRED: The type of blast search to conduct.")
    parser.add_argument("--max_hits", type=int, default=None, help="OPTIONAL: The maximum number of blast matches allowed per input sequence. (May want to set < 300)")
    parser.add_argument("-m", "--merge_strategy", required=False, default="span", choices=["span", "nospan", "all"], help="DEFAULT = 'span'. The strategy for dealing with multiple non-overlapping blast coordinates. See documentation for details.")

    return parser.parse_args()

def delimit(iterable,splitstring):
    return [list(g) for k,g in itertools.groupby(iterable,lambda x:x in splitstring) if not k]

def make_nested_dir(in_dir, string):
    dirname = in_dir+'/'+string
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    return dirname

def find_files_diff_dir(main_dir, search_dir, string):
    os.chdir(search_dir)
    f_list = sorted([f for f in os.listdir('.') if f.endswith(string)])
    os.chdir(main_dir)
    return f_list

def split_name(string, index, delimiter):
    index = int(index)
    name = string.split(delimiter)[index]
    return name

def get_seq_dict(search_dir, end_dir, string):
    os.chdir(search_dir)
    for f in os.listdir('.'):
        if f.startswith(string):
            record_dict = SeqIO.index(f, "fasta")
    os.chdir(end_dir)
    return record_dict
    
def Cluster(fasta, out_dir):
    names = fasta.split('.')
    call_string = "cd-hit-est -i {0}.fasta -o {0}_Out -c 0.8 -n 4 -M 16000".format(names[0])
    proc = sp.call(call_string, shell=True)
    print '\n\n\n'

    for f in os.listdir('.'):
        if f.startswith("{}_Out".format(names[0])):
            shutil.move(f, out_dir)
        
def Parse_Clusters(cluster_dir, cluster_file, out_dir, main_dir):
    os.chdir(cluster_dir)
    #get name of gene from file
    target = cluster_file.split('.')[0]
    lines = []
    print 'Parsing clusters for {}.\n'.format(cluster_file)
    #A clstr file delimits clusters with a >, then numbers all included accession numbers
    #>Cluster 0
    #0	616nt, >JN881132.1... at +/86.20%
    #1	561nt, >KU765220.1... at +/96.08%
    #2	558nt, >KU765300.1... at +/96.59%
    #3	522nt, >JF818216.1... at *
    #4	407nt, >JN881153.1... at +/84.52%
    #Collect the file contents then delimit lists based on presence of 'Cluster' in text
    with open(cluster_file, 'r') as fh_c:
        for line in fh_c:
            line = line.strip()
            if line.startswith('>'):
                lines.append("Cluster")
            else:
                lines.append(line)
    #use delimit function to parse clusters into a list of lists
    sublists = delimit(lines,("Cluster",))
    
    #get seq recs from full fasta file in dictionary format
    gene_name = cluster_file.split('_')[0]
    gene_dict = get_seq_dict(main_dir, cluster_dir, gene_name)
    
    #Now parse through the cd-hit format lines to find accession numbers, write to new fasta
    cluster_number = int(0)
    for cluster in sublists:
        outfile = "{0}_Cluster_{1}.fasta".format(target,cluster_number)
        with open(outfile, 'a') as fh_out:
            #pull accession records from strings looking like this:
            #0	616nt, >JN881132.1... at +/86.20%
            #However sometimes gives odd or incomplete results like this:
            #0	10
            #So test if cluster has actual information before parsing info
            for string in cluster:
                if '>' in string and '...' in string:
                    string1 = string.split('>')
                    string2 = string1[1].split('...')
                    accession = string2[0]
                    #use SeqIO indexed fasta file to pull record and write to
                    #cluter-number-labeled fasta file
                    fh_out.write((gene_dict[accession]).format("fasta"))
        shutil.move(outfile, out_dir)
        cluster_number += 1
        
    #Now grab last sequences in the other outfile from cd-hit-est, write to final cluster file
    #These should always appear in the .clstr file but there may be possible exceptions
    #These should always be duplicates of seqs in the parsed clusters, but this is handled downstream
    records = list(SeqIO.parse(target,'fasta'))
    outfile = "{0}_Cluster_{1}.fasta".format(target,cluster_number)
    with open(outfile, 'a') as fh_out:
        for rec in records:
            fh_out.write(">{}\n{}\n".format(rec.description, str(rec.seq)))
    shutil.move(outfile, out_dir)
                
def Identify_Clusters(in_dir, main_dir):
    os.chdir(in_dir)
    gene_set = set()
    for f in os.listdir('.'):
        gene_set.add(f.split('_')[0])
    gene_list = list(gene_set)
    gene_list.sort()
    
    cluster_info = []
    for gene in gene_list:
        gene_files = []
        for f in os.listdir('.'):
            if f.split('_')[0] == gene:
                with open(f, 'r') as fh_temp:
                    rec_count = int(0)
                    for line in fh_temp:
                        if line.startswith('>'):
                            rec_count += 1
                    gene_files.append([f, rec_count])
        gene_files.sort(key=operator.itemgetter(1), reverse=True)
        cluster_info.append(gene_files)
    os.chdir(main_dir)
    return cluster_info

def make_db(blast_db):
    mdb_str =  "makeblastdb -in {} -dbtype nucl".format(blast_db)
    print '\n\n', mdb_str
    proc = sp.call(mdb_str, shell=True)

def blastn_to_db(task, blast_db, emp_fasta, outname, max_hits):
    if max_hits is None:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 > {3}".format(task, blast_db, emp_fasta, outname)
    	print blast_str
    else:
    	blast_str = "blastn -task {0} -db {1} -query {2} -outfmt 6 -max_target_seqs {3} > {4}".format(task, blast_db, emp_fasta, max_hits, outname)
    	print blast_str
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

def span_coords(coords):
    '''
    Check if multiple bp intervals are within 100bp of each
    other, and if so, merge them and return a new list
    of updated coordinates.
    '''
    new_coords = []
    for i in range(0, (len(coords) - 1)):
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
    print "Parsing contents in {}...\n".format(f)
    accn_set = set()
    file_contents = []
    with open(f, 'r') as fh_in:
        for line in fh_in:
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
            #find accn match, but exclude self-hits
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
    print "\t\tFinished retrieving blast coordinates for {} accessions.\n".format(count)
    
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
    print "\t\tWrote a total of {0} sequences to {1}\n".format(count, autoname) 

    
def Cluster_Blast(in_dir, cluster_info, blast_dir, fasta_list, main_dir, trim_dir, task, max_hits, merge_strategy):
    #move in to parsing directory
    os.chdir(in_dir)
    #iterate over gene sublists
    for gene_group in cluster_info:
        tb = datetime.now()
        #get gene name by splitting first cluster file name
        gene_name = split_name(gene_group[0][0], 0, '_')
        print '\n\n\n--------------------------------------------------------------'
        print "Processing clusters for {}".format(gene_name)
        print '--------------------------------------------------------------'
        
        #make blast database
        #this will come from the largest size cluster file (presumably with actual gene sequences)
        reference = gene_group[0][0]
        make_db(reference)

        #iterate over all gene cluster files and blast to ref
        for cluster_file_list in gene_group:
            outname = cluster_file_list[0].split('.')[0]+"_blast_results.txt"
            blastn_to_db(task, reference, cluster_file_list[0], outname, max_hits)

        #combine blast results into single file for easier writing of new fasta
        merged_blast = "{}_blast_results_merged.txt".format(gene_name)
        with open(merged_blast, 'a') as fh_blast:
            for f in os.listdir('.'):
                if f.startswith(gene_name) and f.endswith("_blast_results.txt"):
                    with open(f, 'r') as fh_temp:
                        for line in fh_temp:
                            fh_blast.write(line)

        #parse blast results
        parsing_list = parse_blastn_output6_NoSelfHits(merged_blast, merge_strategy)

        #cleanup blast results
        for f in os.listdir('.'):
            if f.endswith("_blast_results_merged.txt") or f.endswith("_blast_results.txt") or f.endswith(".nhr") or f.endswith(".nin") or f.endswith(".nsq"):
                shutil.move(f, blast_dir)

        #locate the name of the original gene fasta file
        for f in fasta_list:
            if f.startswith(gene_name):
                emp_fasta = f
        #move to main directory and write sequences (based on coordinates) to fasta file
        os.chdir(main_dir)
        pull_records(emp_fasta, parsing_list)
        
        #move outputs to trimmed directory
        for f in os.listdir('.'):
            if f.endswith("_extracted.fasta") or f.startswith("Log_File_"):
                shutil.move(f, trim_dir)
            
        #move back to parsing directory
        os.chdir(in_dir)

        #show elapsed time
        tf = datetime.now()
        te = tf - tb
        print "\n---------------------------------------------------------------------------\n"
        print "Total time to process {0}: {1} (H:M:S)\n\n".format(gene_name,te)

        
#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    tb = datetime.now()
    os.chdir(args.in_dir)

    cluster_dir = make_nested_dir(args.in_dir, '01_Clustering_Results')
    fasta_list = find_files_diff_dir(args.in_dir, args.in_dir, '.fasta')
    
    print '\n\n\n=====================================================\n'
    print "Performing cd-hit-est to cluster sequences\n"
    print '=====================================================\n\n'
    for f in fasta_list:
        Cluster(f, cluster_dir)

    parse_dir = make_nested_dir(args.in_dir, '02_Parsed_Results')
    cluster_list = find_files_diff_dir(args.in_dir, cluster_dir, '.clstr')
    
    print '\n\n\n======================================================\n'
    print "Parsing results and writing clustered fasta files\n"
    print '======================================================\n\n'
    for f in cluster_list:
        Parse_Clusters(cluster_dir, f, parse_dir, args.in_dir)
    os.chdir(args.in_dir)

    cluster_info = Identify_Clusters(parse_dir, args.in_dir)
    print '\n\n\n====================================================\n'
    print "Summary of clusters generated for all loci\n"
    print '======================================================\n\n'

    for c in cluster_info:
        print "\n{} Clusters:".format(c[0][0].split('_')[0])
        for d in c:
            print "{0} sequences in {1}".format(d[1],d[0])
    
    blast_dir = make_nested_dir(args.in_dir, '03_Blast_Results')
    trim_dir = make_nested_dir(args.in_dir, '04_Trimmed_Results')
    print '\n\n\n====================================================\n'
    print "Performing blasting, merging, and sequence extraction\n"
    print '======================================================\n\n'
    Cluster_Blast(parse_dir, cluster_info, blast_dir, fasta_list, args.in_dir, trim_dir, args.blast_task,args.max_hits,args.merge_strategy)
    
    tf = datetime.now()
    te = tf - tb
    print "\n\nTotal time to completion: {0} (H:M:S)\n\n".format(te)
    
if __name__ == '__main__':
    main()
