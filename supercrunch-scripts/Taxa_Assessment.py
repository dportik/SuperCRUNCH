'''
SuperCRUNCH: Taxa_Assessment module

    Taxa_Assessment: The goal of this script is to examine a fasta file of 
    sequences to see if name labels in the decription line of sequence records can be
    matched those from a taxonomic database. Taxon names can contain a mix of species 
    (two-part name) and subspecies (three-part name) labels. Note that 'subspecies' refers 
    to a three-part name, where the third part can be an actual subspecies label or a 
    unique identifier (such as fied/museum code, or alpha-numerical code). The user supplies 
    a text file containing a list of taxon names to cross-reference. The taxon names in
    the supplied database and the description lines are pre-processed before searches and 
    as a result are NOT case-sensitive. If records with valid 'subspecies' (three part) names 
    should not be included, use the --no_subspecies flag. If this flag is used, only the 
    genus and species will be included from the taxon names database, and only the genus 
    and species will be searched in the sequence records. Output files are written to the 
    output directory specified.
    
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
import time
import sqlite3
import numpy as np
from Bio import SeqIO
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Taxa_Assessment: The goal of this script is to examine a fasta file of 
    sequences to see if name labels in the decription line of sequence records can be
    matched those from a taxonomic database. Taxon names can contain a mix of species 
    (two-part name) and subspecies (three-part name) labels. Note that 'subspecies' refers 
    to a three-part name, where the third part can be an actual subspecies label or a 
    unique identifier (such as fied/museum code, or alpha-numerical code). The user supplies 
    a text file containing a list of taxon names to cross-reference. The taxon names in
    the supplied database and the description lines are pre-processed before searches and 
    as a result are NOT case-sensitive. If records with valid 'subspecies' (three-part) names 
    should not be included, use the --no_subspecies flag. If this flag is used, only the 
    genus and species will be included from the taxon names database, and only the genus 
    and species will be searched in the sequence records. Output files are written to the 
    output directory specified.

    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to a fasta file of "
                            " sequence data.")
    
    parser.add_argument("-t", "--taxa",
                            required=True,
                            help="REQUIRED: The full path to a text file containing "
                            "all taxon names to cross-reference in the fasta file.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to "
                            "write output files.")
    
    parser.add_argument("--no_subspecies",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Ignore any subspecies labels in both the name "
                            "database and record searches (only search binomial names).")
    
    parser.add_argument("--sql_db",
                            required=False,
                            help="OPTIONAL: The full path to the sql database to use "
                            "for searches. Assumes the database was created with "
                            "this module for the input fasta file being used.")
    
    return parser.parse_args()

def parse_taxa(f):
    """
    Retrieve taxon names from user supplied taxon names file (f).
    Will find species and subspecies names, and account for 
    species names that are only represented in subspecies labels.
    Summarizes number of each and will return separate lists 
    for species and subspecies, with all names in uppercase.
    """
    fname = f.split('/')[-1]
    print("\n--------------------------------------------------------------------------------------\n")
    print("Parsing taxon information from: '{}'".format(fname))
    
    with open(f, 'r') as fh:
        species_set = set([line.upper().strip() for line in fh
                               if len(line.split()) == int(2)])
        
    with open(f, 'r') as fh:
        subspecies_set = set([line.upper().strip() for line in fh
                                  if len(line.split()) == int(3)])
        
    #Sometimes a species (binomial) isn't actually present
    #and is only contained in a subspecies (trinomial) name.
    #Just in case, extract the binomial from the trinomial
    #name and add to species set.
    nominates = set([" ".join(x.split()[0:2]) for x in subspecies_set])
    species_set.update(nominates)
    
    #convert sets to sorted lists, return
    species = sorted(species_set)
    subspecies = sorted(subspecies_set)
    
    print("\tFound {:,} species names and {:,} subspecies names."
                  .format(len(species), len(subspecies)))
        
    return species, subspecies


def test_name(name_list, name_parts, test):
    """
    Construct a name from a list containing several
    strings (name_parts), and if option 'test' is desired 
    check if that constructed name is present in a list (name_list).
    """
    if name_parts:
        joined = " ".join(name_parts)
        if test is True:
            if joined in name_list:
                name = joined
            else:
                name = "NA"
        else:
            name = joined
    else:
        name = "NA"

    return name

def parse_fasta_record(line, species, subspecies):
    """
    Deconstruct a record line from a fasta file into
    several elements which will be added to the SQL
    database. Returns an accession number, species label,
    and subspecies label.Name labels are constructed 
    based on several arguments.
    """
    #get accession number from line
    accession = [l.strip('>') for l in line.split()][0]

    #convert line to uppercase for other strings
    line = line.upper()
    
    #get the 'species' name - first two strings following
    #the accession number
    sp_parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                  if len(line.split()) >= int(3)][1:3]
    #construct name and test if it is in the taxon list (species),
    #generates actual name or "NA"
    sp_name = test_name(species, sp_parts, False)
    
    #get the 'subspecies' name - first three strings following
    #the accession number
    ssp_parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                      if len(line.split()) >= int(4)][1:4]
    #construct name and test if it is in the taxon list (subspecies),
    #generates actual name or "NA"
    ssp_name = test_name(subspecies, ssp_parts, False)
        
    return accession, sp_name, ssp_name

def build_sql_db(f, species, subspecies):
    """
    Add the entire contents of this fasta file (f)
    to the SQL database.
    """
    curpath = os.getcwd()
    db = os.path.join(curpath, "Taxa-Assessment.sql.db")
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
        sspname text NOT NULL)""")
    #cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    cur = conn.cursor()
    sql_add_records = "INSERT INTO records VALUES (?, ?, ?)"

    prefix = f.split('.')[0]
    print("\tAdding {} to SQL database...".format(prefix))

    lcnt = int(0)
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                lcnt += 1
                if lcnt % 5000 == 0:
                    print("\tAdded {:,} records...".format(lcnt))
                acc, sp, ssp = parse_fasta_record(line, species, subspecies)
                cur.execute(sql_add_records, (acc, sp, ssp))
                
    f = datetime.now()
    e = f - b
    print("\nTotal time to build SQL database: {0} (H:M:S)\n".format(e))
    
    conn.commit()
    conn.close()

    return db

def run_query(cur, whichname, qseq, qlist):
    
    acc_set, name_set = set(), set()
    
    if whichname == "spname":
        sql_query = "SELECT * FROM records WHERE spname IN ({0})".format(qseq)
        cur.execute(sql_query, qlist)
        results = cur.fetchall()
        [acc_set.add(r['accession']) for r in results]
        [name_set.add(r['spname']) for r in results]
        
    elif whichname == "sspname":
        cur.execute("SELECT * FROM records WHERE sspname IN ({0})'".format(qseq))
        cur.execute(sql_query, qlist)
        results = cur.fetchall()
        [acc_set.add(r['accession']) for r in results]
        [name_set.add(r['sspname']) for r in results]

    return acc_set, name_set


def match_taxa(db, species, subspecies, no_subspecies):
    print("\n--------------------------------------------------------------------------------------\n")
    print("Beginning taxon name matching.")
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    set_accs, setsp, setssp = set(), set(), set()
    
    spseq = ','.join(['?']*len(species))
    sspseq = ','.join(['?']*len(subspecies))

    print("\n\tStarting SQL queries...")
    namecount = int(0)
    if no_subspecies:
        print("\tSearching for species names...")
        accs, spnames = run_query(cur, "spname", spseq, species)
        set_accs.update(accs)
        setsp.update(spnames)
        print("\t\tFinished.")
    
    else:
        print("\tSearching for species names...")
        accs, spnames = run_query(cur, "spname", spseq, species)
        set_accs.update(accs)
        setsp.update(spnames)
        print("\t\tFinished.")

        if subspecies:
            print("\tSearching for subspecies names...")
            accs, sspnames = run_query(cur, "sspname", sspseq, subspecies)
            set_accs.update(accs)
            setssp.update(sspnames)
            print("\t\tFinished.")

    matched_accs = sorted(set_accs)
    matched_spnames = sorted(setsp)
    matched_sspnames = sorted(setssp)

    print("\n\tFound {:,} sequences with matched taxon names.".format(len(matched_accs)))
    print("\tFound {:,} unique matched species (two-part) names.".format(len(matched_spnames)))
    if not no_subspecies:
        print("\tFound {:,} unique matched subspecies (three-part) names.".format(len(matched_sspnames)))
        
    conn.close()

    return matched_accs, matched_spnames, matched_sspnames

def get_unmatched_accs_names(db, matched_accs, no_subspecies):
    print("\n--------------------------------------------------------------------------------------\n")
    print("Gathering records with unmatched names.\n")
    print("\tStarting SQL query...")
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    badaccs, badsp, badssp = set(), set(), set()
    
    accseq = ','.join(['?']*len(matched_accs))
    sql_query = ("SELECT * FROM records WHERE accession NOT IN ({0})".format(accseq))
    cur.execute(sql_query, matched_accs)
    
    results = cur.fetchall()
    [badaccs.add(r['accession']) for r in results]
    [badsp.add(r['spname']) for r in results]
    [badssp.add(r['sspname']) for r in results]
        
    print("\t\tFinished.\n\n\tGathering sequences...")
    
    unmatched_accs = sorted(badaccs)
    unmatched_spnames = sorted(badsp)
    unmatched_sspnames = sorted(badssp)
    
    print("\n\tFound {:,} sequences with an unmatched taxon name.".format(len(unmatched_accs)))
    print("\tFound {:,} unique unmatched species (two-part) names.".format(len(unmatched_spnames)))
    if not no_subspecies:
        print("\tFound {:,} unique unmatched subspecies (three-part) names.".format(len(unmatched_sspnames)))
    print("\n--------------------------------------------------------------------------------------\n")
    
    conn.close()
    
    return unmatched_accs, unmatched_spnames, unmatched_sspnames

def write_log(inlist, outname, namecase):
    if namecase:
        with open(outname, 'a') as fh:
            for x in inlist:
                fh.write("{}\n".format(x.capitalize()))
    else:
        with open(outname, 'a') as fh:
            for x in inlist:
                fh.write("{}\n".format(x))

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
    
def write_fasta(records, acclist, outname):
    tb = datetime.now()
    print("\tWriting {:,} sequence records to {}...".format(len(acclist), outname))
    with open(outname, 'a') as fh:
        for a in acclist:
            fh.write((records[a]).format("fasta"))
    tf = datetime.now()
    te = tf - tb
    print("\t\tElapsed time: {0} (H:M:S)\n\n".format(te))
    
def main():
    args = get_args()
    tb = datetime.now()
    
    os.chdir(args.outdir)
    
    species, subspecies = parse_taxa(args.taxa)
    
    if not args.sql_db:
        db = build_sql_db(args.input, species, subspecies)
    else:
        db = args.sql_db
        print("\n--------------------------------------------------------------------------------------\n")
        print("Using SQL database provided: {}\n".format(db.split('/')[-1]))    
    
    matched_accs, matched_spnames, matched_sspnames =  match_taxa(db, species, subspecies, args.no_subspecies)
    allmatchednames = sorted(matched_spnames + matched_sspnames)
    
    unmatched_accs, unmatched_spnames, unmatched_sspnames = get_unmatched_accs_names(db, matched_accs, args.no_subspecies)
    
    write_log(matched_accs, "Matched_Records_Accession_Numbers.log", False)
    write_log(unmatched_accs, "Unmatched_Records_Accession_Numbers.log", False)
    write_log(allmatchednames, "Matched_Records_Taxon_Names.log", True)
    write_log(unmatched_spnames, "Unmatched_Records_Species_Names.log", True)
    if args.no_subspecies:
        pass
    else:
        write_log(unmatched_sspnames, "Unmatched_Records_Subspecies_Names.log", True)

    records = index_fasta(args.input)
    write_fasta(records, matched_accs, "Matched_Taxa.fasta")
    write_fasta(records, unmatched_accs, "Unmatched_Taxa.fasta")
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    

if __name__ == '__main__':
    main()
