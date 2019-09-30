''''
SuperCRUNCH: Parse_Loci module

    Parse_Loci: This module can be used to find sequences in a fasta file based on a 
    list of taxon names (-t) and a list of locus search terms (-l), and write those sequences
    to locus-specific fasta files. For a sequence to be written to a locus-specific fasta 
    file, it must match either the gene abbreviation or description for that locus AND have 
    a taxon label that is present in the taxon names list. The taxon names list can contain 
    a mix of species (two-part) and subspecies (three-part) names. The 
    --no_subspecies flag can be used to only include species names in searches. In this case, 
    only species names will be considered for records, regardless of whether or not they 
    have a valid subspecies name. 

    All searches occur using SQL, and an SQL database is constructed from the input file (-i) 
    in the output directory specified (-o). If the database for the input file has already been 
    made, the full path to it can be specified using the --sql_db flag, which will save time for
    multiple runs on very large fasta files. All output fasta files and a summary file are 
    written to their relevant directories within the output directory specified (-o). 

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
from Bio import SeqIO
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------    
    Parse_Loci: This module can be used to find sequences in a fasta file based on a 
    list of taxon names (-t) and a list of locus search terms (-l), and write those sequences
    to locus-specific fasta files. For a sequence to be written to a locus-specific fasta 
    file, it must match either the gene abbreviation or description for that locus AND have 
    a taxon label that is present in the taxon names list. The taxon names list can contain 
    a mix of species (two-part) and subspecies (three-part) names. The 
    --no_subspecies flag can be used to only include species names in searches. In this case, 
    only species names will be considered for records, regardless of whether or not they 
    have a valid subspecies name. 

    All searches occur using SQL, and an SQL database is constructed from the input file (-i) 
    in the output directory specified (-o). If the database for the input file has already been 
    made, the full path to it can be specified using the --sql_db flag, which will save time for
    multiple runs on very large fasta files. All output fasta files and a summary file are 
    written to their relevant directories within the output directory specified (-o). 

    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--input",
                            required=True,
                            help="REQUIRED: The full path to a fasta file of "
                            "sequence data.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to "
                            "write output files.")
    
    parser.add_argument("-l", "--loci",
                            required=True,
                            help="REQUIRED: The full path to a text file containing "
                            "loci information to search for within the fasta file.")
    
    parser.add_argument("-t", "--taxa",
                            required=True,
                            help="REQUIRED: The full path to a text file containing all "
                            "taxon names to cross-reference in the fasta file.")
    
    parser.add_argument("--no_subspecies",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Ignore subspecies labels in searches and "
                            "only write species names in the updated description lines "
                            "for sequences.")
    
    parser.add_argument("--exclude",
                            required=False,
                            help="OPTIONAL: The full path to a text file containing a "
                            "list of accession numbers to ignore during searches.")
    
    parser.add_argument("--sql_db",
                            required=False,
                            help="OPTIONAL: The full path to the sql database to use "
                            "for searches. Assumes the database was created with "
                            "this module for the input fasta file being used.")
    
    return parser.parse_args()

def parse_loci_terms(f):
    """
    Input file is tab delimited, three columns:
    [Locus ID] [Locus abbreviation] [Locus string]
    Columns two and three can contain multiple search items,
    but if so must be separated by semicolon character (;).
    """
    fname = f.split('/')[-1]
    print("\n--------------------------------------------------------------------------------------\n")
    print("Parsing locus information from: '{}'".format(fname))
    loci_info = []
    
    with open(f, 'r') as fh_f:
        for line in fh_f:
            #if line not blank, split by tab, replace all quotes and remove whitespace
            cols = [l.replace('\"','').strip().upper() for l in line.split('\t')
                        if line.strip()]
            #for sublist in list, split by semi-colon
            split_cols = [c.split(";") for c in cols]
            loci_info.append(split_cols)
            
    print("\tFound {:,} loci to search.".format(len(loci_info)))
    return loci_info

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

def index_fasta(f):
    """
    Use SeqIO to index the fasta file (f), which may be too big 
    to parse into a list. Returns a dictionary-like structure.
    """
    b = datetime.now()
    fname = f.split('/')[-1]
    print("\n--------------------------------------------------------------------------------------\n")
    print("Indexing fasta file: '{}'\n\tThis could take some time...".format(fname))
    
    records = SeqIO.index(f, "fasta")
    
    f = datetime.now()
    e = f - b
    print("\nTotal time to index fasta file: {0} (H:M:S)\n".format(e))

    return records

def test_name(name_list, name_parts, test):
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

def parse_fasta_record(line, species, subspecies, no_subspecies):
    #get accession number from line
    accession = [l.strip('>') for l in line.split()][0]

    #convert line to uppercase for other strings
    line = line.upper()
    
    #get the 'species' name - first two strings following
    #the accession number
    sp_parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                  if len(line.split()) >= int(3)][1:3]
    #construct name but do not test if in taxon list (species);
    #produces binomial name or NA
    sp_name = test_name(species, sp_parts, False)

    #decide whether to attempt subspecies name construction
    #if no subspecies desired, make all subspecies labels "NA"
    if no_subspecies is True:
        ssp_name = "NA"

    #if subspecies desired
    else:
        #get the 'subspecies' name - first three strings following
        #the accession number
        ssp_parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                      if len(line.split()) >= int(4)][1:4]
        #construct name and test if it is in the taxon list (subspecies),
        #generates actual name or "NA"
        ssp_name = test_name(subspecies, ssp_parts, True)

    #get the 'description' line, which here is everything after the 'species' name
    #check to make sure there was actually a species name and stuff follows
    if len(line.split()[3:]) >= int(1):
        description = (" ".join(line.replace(",", '')
                         .replace(";", '')
                         .replace(":", '')
                         .replace(")", '')
                         .replace("(", '')
                         .replace("<", '')
                         .split()[3:]))+" "
    else:
        description = "NA"

    #try to obtain a field/museum/sample code using the keywords
    #voucher, isolate, and strain in the description line
    if 'VOUCHER' in description:
        parts = description.split('VOUCHER ')[-1].split()
        if parts[0].replace("-", "").isalpha() and len(parts) > 1:
            voucher = "Voucher_{}_{}".format(parts[0], parts[1])
        else:
            voucher = "Voucher_{}".format(parts[0])
            
    elif 'ISOLATE' in description:
        parts = description.split('ISOLATE ')[-1].split()
        if parts[0].replace("-", "").isalpha() and len(parts) > 1:
            voucher = "Voucher_{}_{}".format(parts[0], parts[1])
        else:
            voucher = "Voucher_{}".format(parts[0])
            
    elif 'STRAIN' in description:
        parts = description.split('STRAIN ')[-1].split()
        if parts[0].replace("-", "").isalpha() and len(parts) > 1:
            voucher = "Voucher_{}_{}".format(parts[0], parts[1])
        else:
            voucher = "Voucher_{}".format(parts[0])
    else:
        voucher = "NA"
        
    return accession, sp_name, ssp_name, description, voucher

def build_sql_db(f, db, species, subspecies, no_subspecies):
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
        description text NOT NULL,
        voucher text NOT NULL)""")
    #cur.execute("SELECT name FROM sqlite_master WHERE type='table'")

    sql_add_records = "INSERT INTO records VALUES (?, ?, ?, ?, ?)"

    lcnt = 0
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                lcnt += 1
                if lcnt % 5000 == 0:
                    print("\tAdded {:,} records...".format(lcnt))
                acc, sp, ssp, d, v = parse_fasta_record(line, species, subspecies, no_subspecies)
                cur.execute(sql_add_records, (acc, sp, ssp, d, v))
    conn.commit()
    
    #cur.execute("SELECT accession, spname, sspname, description, voucher FROM records")
    #print cur.fetchall()
    conn.close()
    f = datetime.now()
    e = f - b
    print("\nTotal time to build SQL database: {0} (H:M:S)\n".format(e))


def prep_terms(l):
    padded = [" {} ".format(t) for t in l[1]]
    terms = padded + l[2]
    return terms

def perform_searches(cur, l, spseq, species, exclude_accs):
    print("\nSearching for {}:\n".format(l[0][0]))
    
    terms = prep_terms(l)

    acc_set = set()
    excluded = int(0)
    
    for t in terms:
        sql_query = "SELECT * FROM records WHERE spname IN ({0}) AND description LIKE '%{1}%'".format(spseq, t)
        #print sql_query
        cur.execute(sql_query, species)
        results = cur.fetchall()
        print("\t{:,} records found using term: '{}'".format(len(results), t))
        for r in results:
            if r['accession'] not in exclude_accs:
                acc_set.add(r['accession'])
            else:
                excluded += 1
    if excluded >= 1:
        print("\n\t*Excluded {} records for {}, based on --exclude flag.".format(excluded, l[0][0]))
    print("\n\t{:,} total unique records found for {}.".format(len(acc_set), l[0][0]))


    return sorted(acc_set)
    
def get_records(gene, cur, accessions, records):
    accseq = ','.join(['?']*len(accessions))
    sql_query = "SELECT * FROM records WHERE accession IN ({})".format(accseq)
    cur.execute(sql_query, accessions)
    results = cur.fetchall()

    return results

def write_fasta(gene, writing_recs, records, no_subspecies):
    #accession, spname, sspname, description, voucher
    outname = "{}.fasta".format(gene)
    print("\tWriting {:,} records to {}".format(len(writing_recs), outname))
    
    if no_subspecies is True:
        with open(outname, 'a') as fh:
            for r in writing_recs:
                if r['voucher'] != 'NA':
                    fh.write(">{0} {1} {2} DESCRIPTION {3}\n{4}\n"
                                 .format(r['accession'],
                                             r['spname'].capitalize(),
                                             r['voucher'],
                                             r['description'].lower(),
                                             records[r['accession']].seq))
                elif r['voucher'] == 'NA':
                    fh.write(">{0} {1} DESCRIPTION {2} \n{3}\n"
                                 .format(r['accession'],
                                             r['spname'].capitalize(),
                                             r['description'].lower(),
                                             records[r['accession']].seq))
                    
    elif no_subspecies is False:
        with open(outname, 'a') as fh:
            for r in writing_recs:
                if r['sspname'] == 'NA':
                    if r['voucher'] != 'NA':
                        fh.write(">{0} {1} {2} DESCRIPTION {3}\n{4}\n"
                                     .format(r['accession'],
                                                 r['spname'].capitalize(),
                                                 r['voucher'],
                                                 r['description'].lower(),
                                                 records[r['accession']].seq))
                    elif r['voucher'] == 'NA':
                        fh.write(">{0} {1} DESCRIPTION {2} \n{3}\n"
                                     .format(r['accession'],
                                                 r['spname'].capitalize(),
                                                 r['description'].lower(),
                                                 records[r['accession']].seq))
        
                elif r['sspname'] != 'NA':
                    if r['voucher'] != 'NA':
                        fh.write(">{0} {1} {2} DESCRIPTION {3}\n{4}\n"
                                     .format(r['accession'],
                                                 r['sspname'].capitalize(),
                                                 r['voucher'],
                                                 r['description'].lower(),
                                                 records[r['accession']].seq))
                    elif r['voucher'] == 'NA':
                        fh.write(">{0} {1} DESCRIPTION {2} \n{3}\n"
                                     .format(r['accession'],
                                                 r['sspname'].capitalize(),
                                                 r['description'].lower(),
                                                 records[r['accession']].seq))
    
def write_log(loci_counts):
    outname = 'Loci_Record_Counts.txt'
    with open(outname, 'a') as fh:
        fh.write("Locus_Name\tRecords_Written\n")
        for l in loci_counts:
            fh.write("{}\t{}\n".format(l[0], l[1]))
    
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
    fdir = os.path.join(curpath, "Parsed-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)
        
    sdir = os.path.join(curpath, "Summary-File")
    if not os.path.exists(sdir):
        os.mkdir(sdir)
        

    return fdir, sdir

def cleanup(fdir, sdir):
    """
    Moves relevant output files to their output directories. 
    """
    print("\n\nCleaning up output files...")
    
    [os.remove(f) for f in os.listdir('.') if f.endswith('.fasta') and os.stat(f).st_size == 0]
    [shutil.move(f, fdir) for f in os.listdir('.') if f.endswith('.fasta')]
    shutil.move('Loci_Record_Counts.txt', sdir)
    print("\tDone!")
        
def search_runner(species, loci_info, db, records, no_subspecies, exclude):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    
    spseq = ','.join(['?']*len(species))

    loci_counts = []

    if exclude:
        with open(exclude, 'r') as fh:
            exclude_accs = [l.strip() for l in fh]
    else:
        exclude_accs = []

    print("\n--------------------------------------------------------------------------------------\n")
    for l in loci_info:
        b = datetime.now()
        
        accessions = perform_searches(cur, l, spseq, species, exclude_accs)
        writing_recs = get_records(l[0][0], cur, accessions, records)
        write_fasta(l[0][0], writing_recs, records, no_subspecies)
        loci_counts.append([l[0][0], len(writing_recs)])

        f = datetime.now()
        e = f - b
        print("\n\tElapsed time: {} (H:M:S).\n".format(e))
    
    write_log(loci_counts)

def main():
    tb = datetime.now()
    args = get_args()
    
    loci_info = parse_loci_terms(args.loci)
    species, subspecies = parse_taxa(args.taxa)

    os.chdir(args.outdir)
    fdir, sdir = make_dirs(args.outdir)
    
    if not args.sql_db:
        prefix = args.input.split('/')[-1].split('.')[0]
        curpath = os.getcwd()
        sql_db_name = os.path.join(curpath, "{}.sqlite.db".format(prefix))
        build_sql_db(args.input, sql_db_name, species, subspecies, args.no_subspecies)
        
    else:
        sql_db_name = args.sql_db
        print("\n--------------------------------------------------------------------------------------\n")
        print("Using SQL database provided: {}\n".format(sql_db_name.split('/')[-1]))
        
    records = index_fasta(args.input)

    search_runner(species, loci_info, sql_db_name, records, args.no_subspecies, args.exclude)
    cleanup(fdir, sdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal time to parse all loci: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")
    
    
if __name__ == '__main__':
    main()
