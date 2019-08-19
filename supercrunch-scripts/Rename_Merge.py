'''
SuperCRUNCH: Rename_Merge module
                            
    Rename_Merge: This module can be used to relabel taxon names that did not match 
    the taxon list used during the Taxa_Assessment step. These records will have 
    been written to a file called Unmatched_Taxa.fasta, which should be used as
    the input file (-i). Currently, two-part names (e.g., Genus species) and 
    three-part names (e.g., Genus species subspecies) can be replaced using a substitute 
    name of any length. That is, a species name can be replaced using a different 
    species name or a subspecies name. The same is true for replacing subspecies 
    labels. The replacement names file (specified by -r) is a tab-delimited text 
    file with two columns. The first column contains the unmatched name to replace, 
    and the second column contains the replacement name. All successfully relabeled 
    records are written to a fasta file called Relabeled.fasta. These records can 
    also be joined with records from an additional fasta file using the -m flag. 
    This is ideal for joining updated records with those from Matched_Taxa.fasta, 
    and will produce an output fasta file called Merged.fasta. Finally, a summary 
    of the number of records relabeled for each name provided is written as 
    Renaming_Summary.txt. All output files are written to the output directory 
    specified (-o). 
    
    
-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
July 2019
Distributed under the 
GNU General Public Lincense
'''
import sys
import os
import argparse
import shutil
import sqlite3
import numpy as np
from collections import Counter
from Bio import SeqIO
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Rename_Merge: This module can be used to relabel taxon names that did not match 
    the taxon list used during the Taxa_Assessment step. These records will have 
    been written to a file called Unmatched_Taxa.fasta, which should be used as
    the input file (-i). Currently, two-part names (e.g., Genus species) and 
    three-part names (e.g., Genus species subspecies) can be replaced using a substitute 
    name of any length. That is, a species name can be replaced using a different 
    species name or a subspecies name. The same is true for replacing subspecies 
    labels. The replacement names file (specified by -r) is a tab-delimited text 
    file with two columns. The first column contains the unmatched name to replace, 
    and the second column contains the replacement name. All successfully relabeled 
    records are written to a fasta file called Relabeled.fasta. These records can 
    also be joined with records from an additional fasta file using the -m flag. 
    This is ideal for joining updated records with those from Matched_Taxa.fasta, 
    and will produce an output fasta file called Merged.fasta. Finally, a summary 
    of the number of records relabeled for each name provided is written as 
    Renaming_Summary.txt. All output files are written to the output directory 
    specified (-o). 
    
    DEPENDENCIES: Python: BioPython, sqlite3.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to a fasta file to replace "
                            "names inside. If you have used the Taxa_Assessment module, "
                            "this should be the file Unmatched_Taxa.fasta.")
    
    parser.add_argument("-r", "--replace",
                            required=True,
                            help="REQUIRED: The full path to a text file containing "
                            "the replacement name information.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-m", "--merge",
                            required=False,
                            default=None,
                            help="OPTIONAL: The full path to a fasta file to merge "
                            "with the updated records. If you have used the "
                            "Taxa_Assessment module, this should be the file "
                            "Matched_Taxa.fasta.")

    parser.add_argument("--sql_db",
                            required=False,
                            help="OPTIONAL: The full path to the sql database to use "
                            "for searches. Assumes the database was created with "
                            "this module for the input fasta file being used.")
    
    parser.add_argument("--quiet",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Show less output while running.")
    
    return parser.parse_args()

def get_relabeling_dict(f):
    """
    Function to parse relabeling file into two distinct
    dictionaries based on whether the name to replace 
    contains two parts (Genus species) or three parts
    (Genus species subspecies). Name to replace is the
    key and the replacement name is the value. Returns 
    both dictionaries.
    """
    fname = f.split('/')[-1]
    print("\n--------------------------------------------------------------------------------------\n")
    print("Parsing contents in: '{}'".format(fname))
    
    #intiate empty dictionaries
    sp_dict, ssp_dict = {}, {}

    #list comprehension to get file contents
    with open(f, 'r') as fh:
        contents = [l.strip().split('\t') for l in fh]

    #test whether name is two or three parts,
    #add to relevant dictionary
    for c in contents:
        if len(c[0].split()) == 3:
            ssp_dict[c[0]] = c[1]
            #print("{} = {}".format(c[0], c[1]))
            
        elif len(c[0].split()) == 2:
            sp_dict[c[0]] = c[1]
            #print("{} = {}".format(c[0], c[1]))
            
    print("\n\tFound {:,} species (two-part) name(s) to replace.".format(len(sp_dict)))
    print("\tFound {:,} subspecies (three-part) name(s) to replace.\n".format(len(ssp_dict)))

    return sp_dict, ssp_dict


def parse_fasta_record(line):
    """
    Deconstruct a record line from a fasta file into
    several elements which will be added to the SQL
    database. Returns an accession number, species label,
    subspecies label, and variations of the remaining 
    description line. The taxon labels are in uppercase.
    """
    #get accession number from line
    accession = [l.strip('>') for l in line.split()][0]

    #get description line assuming a two-part name is present
    twod = " ".join([l for l in line.split() if len(line.split()) >= int(3)][3:])
    if not twod:
        twod = "NA"

    #get description line assuming a three-part name is present
    threed = " ".join([l for l in line.split() if len(line.split()) >= int(4)][4:])
    if not threed:
        threed = "NA"
    
    #convert line to uppercase for other strings
    line = line.upper()
    
    #get the 'species' name - first two strings following
    #the accession number
    sp_name = " ".join([l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                  if len(line.split()) >= int(3)][1:3])
    if not sp_name:
        sp_name = "NA"
    
    #get the 'subspecies' name - first three strings following
    #the accession number
    ssp_name = " ".join([l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                      if len(line.split()) >= int(4)][1:4])
    if not ssp_name:
        ssp_name = "NA"

    #print("{}\n\t{}\n\t{}\n\t{}\n\t{}".format(accession, sp_name, ssp_name, twod, threed))
    return accession, sp_name, ssp_name, twod, threed

def build_sql_db(f, quiet):
    """
    Add the entire contents of this fasta file (f)
    to the SQL database.
    """
    curpath = os.getcwd()
    db = os.path.join(curpath, "Rename-Merge.sql.db")
    b = datetime.now()
    print("\n--------------------------------------------------------------------------------------\n")
    print("Building SQL database: {}".format(db))
    if os.path.exists(db):
        os.remove(db)
        
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    
    cur.execute("DROP TABLE IF EXISTS records;")
    cur.execute("""CREATE TABLE records (
        accession text NOT NULL,
        spname text NOT NULL,
        sspname text NOT NULL,
        twodescription text NOT NULL,
        threedescription text NOT NULL)""")
    #cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    cur = conn.cursor()
    sql_add_records = "INSERT INTO records VALUES (?, ?, ?, ?, ?)"

    lcnt = int(0)
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                lcnt += 1
                if not quiet:
                    if lcnt % 5000 == 0:
                        print("\tAdded {:,} records...".format(lcnt))
                acc, sp, ssp, twod, threed = parse_fasta_record(line)
                cur.execute(sql_add_records, (acc, sp, ssp, twod, threed))
                
    f = datetime.now()
    e = f - b
    print("\nTotal time to build SQL database: {0} (H:M:S)\n".format(e))
    
    conn.commit()
    conn.close()

    return db

def index_fasta(f):
    """
    Use SeqIO index function to load fasta file, which
    is the best option for massive files.
    """
    tb = datetime.now()
    print("\nIndexing fasta file: {}".format(f.split('/')[-1]))
    print("\tThis could take some time...")
    records = SeqIO.index(f, "fasta")
    tf = datetime.now()
    te = tf - tb
    print("\tTotal time to index fasta file: {0} (H:M:S)\n\n".format(te))

    return records

def replace_and_write(sp_dict, ssp_dict, db, fasta, outname, quiet):
    """
    Use sql database (db) to search for names (keys) present in
    ssp_dict and sp_dict in the sequence records. If the name is
    found, a custom description line is created which includes the
    accession number, replacement name, a flag called *Relabeled*, 
    and the remaining original description line (everything after 
    the original taxon label). The sequence is written using the 
    accession number and the biopython dictionary created from
    the fasta file. A list (rename_info) is populated with sublists
    which contain: [original name, replacement name, records replaced],
    and returned. 
    """
    print("\n--------------------------------------------------------------------------------------\n")
    print("Beginning name replacement.")
    b = datetime.now()
    
    #connect to database
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    #build dictionary structure from fasta file
    records = index_fasta(fasta)

    #initiate empty list to store sublists:
    #[original name, replacement name, records replaced]
    rename_info = []
    print("\nSearching names and writing relevant records to: {}.".format(outname))

    #intitiate empty set for accs
    acc_set = set()
    
    #Perform search and replace for three-part names
    #first obtain sorted list of uppercase keys (e.g. names) from ssp_dict
    ssp_keys = sorted([k.upper() for k in ssp_dict.keys()])
    if ssp_keys:
        #create placeholder for sql query, we want to search in list of all ssp names
        sspseq = ','.join(['?']*len(ssp_keys))
        #create query to pull all entries with spname present in ssp_keys, execute
        sql_query = "SELECT * FROM records WHERE sspname IN ({})".format(sspseq)
        cur.execute(sql_query, ssp_keys)
        results = cur.fetchall()
        
        #get all the spname columns from matches, convert to capitalized
        all_ssp_name_matches = [r['sspname'].capitalize() for r in results]
        #count the total number of occurrences of each sspname in the above list
        #using the Counter function, which returns a dict of item:count entries
        #for all items in the list
        ssp_occurences = Counter(all_ssp_name_matches)
        
        #iterate over key val pairs in ssp_occurences dict
        for k in sorted(ssp_occurences.keys()):
            #append [original name, replacement name, number of records replaced] to rename_info
            rename_info.append([k, ssp_dict[k], ssp_occurences[k]])
            if not quiet:
                print("\tRelabeling {:,} records matching name '{}' with updated name '{}'."
                          .format(ssp_occurences[k], k, ssp_dict[k]))
                
        #add accessions to set so we don't also write them for the species match later
        [acc_set.add(r['accession']) for r in results]

        #create list to use for writing updated records from sql query results
        resultlist = [[r['accession'], ssp_dict[r['sspname'].capitalize()], r['threedescription']] for r in results]
        if resultlist:
            #write to output fasta file
            with open(outname, 'a') as fh:
                for r in resultlist:
                    #construct the new description line with accession, replacement name, relabel
                    #flag, remaining description line, and use biopython to get accession in dict
                    #structure to write the sequence
                    fh.write(">{0} {1} {2} {3}\n{4}\n"
                                 .format(r[0], r[1], "*Relabeled*", r[2], records[r[0]].seq))
                    #print(">{0} {1} {2} {3}".format(r[0], r[1], "*Relabeled*", r[2]))

    #Perform search and replace for two-part names
    #first obtain sorted list of uppercase keys (e.g. names) from sp_dict
    sp_keys = sorted([k.upper() for k in sp_dict.keys()])
    if sp_keys:
        #create placeholder for sql query, we want to search in list of all sp names
        spseq = ','.join(['?']*len(sp_keys))
        #create query to pull all entries with spname present in sp_keys, execute
        sql_query = "SELECT * FROM records WHERE spname IN ({})".format(spseq)
        cur.execute(sql_query, sp_keys)
        results = cur.fetchall()

        #get all the spname columns from matches, convert to capitalized
        all_sp_name_matches = [r['spname'].capitalize() for r in results]
        #count the total number of occurrences of each spname in the above list
        #using the Counter function, which returns a dict of item:count entries
        #for all items in the list
        sp_occurences = Counter(all_sp_name_matches)
        
        #iterate over sorted keys in sp_occurences dict
        for k in sorted(sp_occurences.keys()):
            #append [original name, replacement name, number of records replaced] to rename_info
            rename_info.append([k, sp_dict[k], sp_occurences[k]])
            if not quiet:
                print("\tRelabeling {:,} records matching name '{}' with updated name '{}'."
                          .format(sp_occurences[k], k, sp_dict[k]))

        #create list to use for writing updated records from sql query results, but apply additional filter
        #such that no records are written that have already been written using the ssp search
        resultlist = [[r['accession'], sp_dict[r['spname'].capitalize()], r['twodescription']] for r in results
                          if r['accession'] not in acc_set]
            
        if resultlist:
            #write to output fasta file
            with open(outname, 'a') as fh:
                for r in resultlist:
                    #construct the new description line with accession, replacement name, relabel
                    #flag, remaining description line, and use biopython to get accession in dict
                    #structure to write the sequence
                    fh.write(">{0} {1} {2} {3}\n{4}\n"
                                 .format(r[0], r[1], "*Relabeled*", r[2], records[r[0]].seq))
                    #print(">{0} {1} {2} {3}".format(r[0], r[1], "*Relabeled*", r[2]))
    

    #get total number of records replaced by summing third element
    #in rename_info using list comprehension and sum function
    relabeled_count = sum([i[2] for i in rename_info])
    
    f = datetime.now()
    e = f - b
    print("\nRelabeled {0:,} of {1:,} records.".format(relabeled_count, len(records)))
    print("\nTotal time to write relabeled fasta file: {0} (H:M:S)\n".format(e))

    #close db connection
    conn.close()

    return rename_info
            
def merge_fastas(f1, f2, outname):
    """
    Simply add the contents of each fasta file
    (f1, f2) to a combined output file (outname).
    """
    print("\n--------------------------------------------------------------------------------------\n")
    print("Merging updated records with {}.".format(f2.split('/')[-1]))
    b = datetime.now()
    
    with open(outname, 'a') as fhout:
        with open(f1, 'r') as fh1:
            for line in fh1:
                fhout.write(line)
                
        with open(f2, 'r') as fh2:
            for line in fh2:
                fhout.write(line)
                
    f = datetime.now()
    e = f - b
    print("Total time to write merged fasta file: {0} (H:M:S)\n".format(te))
    
def write_log(rename_info):
    """
    Write information in rename_info to an output
    file. Contents are lists: 
    [original name, replacement name, number of records replaced]
    """
    outname = "Renaming_Summary.txt"
    print("\n--------------------------------------------------------------------------------------\n")
    print("Writing summary of relabeling to: {}.\n".format(outname))
    with open(outname, 'a') as fh:
        fh.write("{}\t{}\t{}\n".format("Orig_Name", "Replace_Name", "Records_Relabeled"))
        for r in rename_info:
            fh.write("{}\t{}\t{}\n".format(r[0], r[1], r[2]))
            
def main():
    args = get_args()
    tb = datetime.now()

    sp_dict, ssp_dict = get_relabeling_dict(args.replace)
    
    os.chdir(args.outdir)
    
    if not args.sql_db:
        db = build_sql_db(args.input, args.quiet)
    else:
        db = args.sql_db
        print("\n--------------------------------------------------------------------------------------\n")
        print("Using SQL database provided: {}\n".format(db.split('/')[-1]))
    
    relabeled_fasta, merged_fasta = "Renamed.fasta", "Merged.fasta"
    rename_info = replace_and_write(sp_dict, ssp_dict, db, args.input, relabeled_fasta, args.quiet)
    write_log(rename_info)
    
    if args.merge:
        merge_fastas(relabeled_fasta, args.merge, merged_fasta)
    
    tf = datetime.now()
    te = tf - tb
    print("\n--------------------------------------------------------------------------------------")
    print("\nFinished. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
                 
if __name__ == '__main__':
    main()
