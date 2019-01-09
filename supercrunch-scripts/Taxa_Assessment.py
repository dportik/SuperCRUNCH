'''
SuperCRUNCH: Taxa_Assessment module

Usage: python Taxa_Assessment.py -i [fasta file] REQUIRED
							-t [text file with taxon names to cross-reference] REQUIRED
							-o [path to output directory] REQUIRED
                            --no_subspecies (OPTIONAL flag, excludes subspecies names from searches)
                            
    Taxa_Assessment: The goal of this script is to examine a large fasta file of sequences downloaded from
    GenBank and examine taxon names in the decription line of sequence records to see if
    they match those from a taxonomic database. Taxon names can contain a mix of species (binomial name) 
    and subspecies (trinomial name) labels. The user supplies a text file containing a
    list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should not be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record.
    Two output fasta files are produced, one containing only
    records with taxon names present in the database, and one containing only records which
    contain non-matched taxon names. These fasta files are written to the output directory
    specified. A log file summarizing the matched and unmatched taxon names is also written
    to the output directory, and are called 'Matched_Taxon_Names.log' and
    'Unmatched_Taxon_Names.log'. The unmatched file can form the basis of the file needed
    to rename or correct the taxon ID of these entries in subsequent scripts.

[-t] The input taxon name file should simply contain a list of taxon names, one each line. 
    For all names the genus and species (and subspecies, if present) should separated by a space.
    All taxon names are converted to uppercase for searching, so names are NOT case-sensitive.
 
	Example of file structure for taxon information (showing species and subspecies examples):
	
    Varanus acanthurus
    Varanus albigularis albigularis
    Varanus albigularis microstictus
    Varanus auffenbergi
    Varanus bangonorum
    Varanus baritji
    Varanus beccarii

[-o] Full path to an existing directory in which the outputs will be written. 

	Outputs Files:
	
	Matched_Taxon_Names.log - Text file summarizing each taxon name matched to database. 

    Unmatched_Taxon_Names.log - Text file summarizing each taxon name not matched to database. 
	
	Matched_Taxa.fasta - A fasta file containing all records with a matched taxon name.

	Unmatched_Taxa.fasta - A fasta file containing all records without a matched taxon name.
                           The taxon names in this file can be renamed/corrected in the
                           subsequent script and merged with the Matched_Taxa.fasta. 
                           
[--no_subspecies] Excludes all subspecies names from searches, regardless of whether the
     taxon names file contains them or not. This essentially reduces all subspecies names
     to a binomial name, so all subspecies would be considered a single species. 
    
 -------------------------
For Python 2.7
Python modules required:
	-BioPython (using SeqIO module)
	-NumPy 
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
import argparse
import time
import numpy as np
from Bio import SeqIO
from datetime import datetime

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Taxa_Assessment: The goal of this script is to examine a large fasta file of sequences downloaded from
    GenBank and examine taxon names in the decription line of sequence records to see if
    they match those from a taxonomic database. Taxon names can contain a mix of species (binomial name) 
    and subspecies (trinomial name) labels. The user supplies a text file containing a
    list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should not be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record.
    Two output fasta files are produced, one containing only
    records with taxon names present in the database, and one containing only records which
    contain non-matched taxon names. These fasta files are written to the output directory
    specified. A log file summarizing the matched and unmatched taxon names is also written
    to the output directory, and are called 'Matched_Taxon_Names.log' and
    'Unmatched_Taxon_Names.log'. The unmatched file can form the basis of the file needed
    to rename or correct the taxon ID of these entries in subsequent scripts.

    [-t] The input taxon name file should simply contain a list of taxon names, one each line. 
    For all names the genus and species (and subspecies, if present) should separated by a space.
    All taxon names are converted to uppercase for searching, so names are NOT case-sensitive.

    DEPENDENCIES: Python: BioPython.
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input", required=True, help="REQUIRED: The full path to a fasta file of GenBank sequence data")
    parser.add_argument("-t", "--taxa", required=True, help="REQUIRED: The full path to a text file containing all taxon names to cross-reference in the fasta file.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: The full path to an existing directory to write output files.")
    parser.add_argument("--no_subspecies", required=False, action='store_true', help="OPTIONAL: Ignore any subspecies labels in both the name database and record searches (only search binomial names).")
    return parser.parse_args()

def index_fasta(f):
    tb = datetime.now()
    print "\nIndexing fasta file, this could take some time..."
    record_index = SeqIO.index(f, "fasta")
    tf = datetime.now()
    te = tf - tb
    print "\tTotal time to index fasta file: {0} (H:M:S)\n\n".format(te)

    return record_index

def parse_taxa(f, no_subspecies):
    print "\nParsing taxon information from {}.".format(f)
    with open(f, 'r') as fh_f:
        species_set = set([line.upper().strip() for line in fh_f if len(line.split()) == int(2)])
    with open(f, 'r') as fh_f:
        subspecies_set = set([line.upper().strip() for line in fh_f if len(line.split()) == int(3)])
    #sometimes a species (binomial) isn't present alone and is only in subspecies (trinomial) names
    #extract the binomial from the trinomials and add to species set
    nominate_form = set([" ".join(x.split()[0:2]) for x in subspecies_set])
    species_set.update(nominate_form)
    #convert sets to sorted lists, return
    species = sorted(species_set)
    subspecies = sorted(subspecies_set)
    if no_subspecies is False:
        print "\tFound {} species names and {} subspecies names.\n".format(len(species), len(subspecies))
    elif no_subspecies is True:
        print "\tFound {} species names.\n".format(len(species))
    return species, subspecies

def taxon_match_sp(taxon_sp,species):
    match = False
    if taxon_sp in species:
        match = True
    return match

def taxon_match_ssp(taxon_ssp,subspecies):
    match = False
    if taxon_ssp in subspecies:
        match = True
    return match

def get_taxon(line):
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if line.split() >= int(3)][1:3]
    taxon_sp = " ".join(parts1)
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if line.split() >= int(4)][1:4]
    taxon_ssp = " ".join(parts2)
    return taxon_sp, taxon_ssp

def get_accession(line):
    acc = [l.strip('>') for l in line.split()][0]
    return acc

def search_fasta(f, species, subspecies, no_subspecies):
    acc_match = set()
    acc_unmatch = set()
    tax_match = set()
    tax_unmatch = set()
    
    tb = datetime.now()
    with open(f, 'r') as fh_fasta:
        print "Search progress:"
        lcnt = int(0)
        for line in fh_fasta:
            lcnt += 1
            if lcnt % 100000 == 0:
                print "\tLine {}...".format(lcnt)
            if line.startswith('>'):
                if '|' in line:
                    line.replace("|"," ")
                line = line.upper().strip()
                taxon_sp, taxon_ssp = get_taxon(line)
                acc = get_accession(line)
                #Statement for testing if species and subspecies names
                #occur in the line. This should handle all cases correctly.
                #First set of statements are when subspecies are desired:
                if no_subspecies is False:
                    if taxon_match_sp(taxon_sp,species) is True:
                        if taxon_match_ssp(taxon_ssp,subspecies) is True:
                            acc_match.add(acc)
                            tax_match.add(taxon_ssp.capitalize())
                            #print "{} found in species and {} found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize()) 
                        else:
                            acc_match.add(acc)
                            tax_match.add(taxon_sp.capitalize())
                            #print "{} found in species but {} NOT found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize()) 
                    elif taxon_match_sp(taxon_sp,species) is False:
                        if taxon_match_ssp(taxon_ssp,subspecies) is True:
                            acc_match.add(acc)
                            tax_match.add(taxon_ssp.capitalize())
                            #print "{} NOT found in species but {} found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize())
                        else:
                            acc_unmatch.add(acc)
                            tax_unmatch.add(taxon_sp.capitalize())
                            #print "{} NOT found in species and {} NOT found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize())
                #Second set of statements are when subspecies are desired
                if no_subspecies is True:
                    if taxon_match_sp(taxon_sp,species) is True:
                        acc_match.add(acc)
                        tax_match.add(taxon_sp.capitalize())
                        #print "{} found in species".format(taxon_sp.capitalize()) 
                    elif taxon_match_sp(taxon_sp,species) is False:
                        acc_unmatch.add(acc)
                        tax_unmatch.add(taxon_sp.capitalize())
                        #print "{} NOT found in species".format(taxon_sp.capitalize())
    tf = datetime.now()
    te = tf - tb
    print "Total time to search taxon names in fasta file: {0} (H:M:S)\n\n".format(te)
                
    return acc_match, acc_unmatch, tax_match, tax_unmatch
                        
def write_log(in_set,fname):
    print "Writing results to {}\n".format(fname)
    slist = sorted(in_set)
    with open(fname, 'a') as fh_out:
        for i in slist:
            fh_out.write("{}\n".format(i))

def read_accs(f):
    new_set = set()
    with open(f, 'r') as fh_in:
        for line in fh_in:
            if "|" not in line:
                new_set.add(line.strip())
    return new_set
        
def write_fasta(in_set,prefix,record_index):
    fasta = "{}.fasta".format(prefix)
    print "Writing sequence records to {}".format(fasta)
    tb = datetime.now()
    slist = sorted(in_set)
    upper = find_upper(len(slist))
    acnt = int(0)
    with open(fasta, 'a') as fh_out:
        for a in slist:
            acnt += 1
            interval = upper // 50
            if acnt % interval == 0:
                print "\t{} records written...".format(acnt)
            #reduce line counts with this output format
            fh_out.write(">{0}\n{1}\n".format(record_index[a].description,record_index[a].seq))
            #vs this format from SeqIO
            #fh_out.write((record_index[a]).format("fasta"))
    print "Wrote {0} records to {1}".format(len(in_set),fasta)
    tf = datetime.now()
    te = tf - tb
    print "\tTotal time to write fasta file: {0} (H:M:S)\n\n".format(te)
    
def find_upper(count):
    vala = np.array([100,1000,10000])
    valb = np.arange(100000,1000000000,100000)
    vals = np.append(vala,valb)
    for i in vals:
        if count <= i:
            upper = i
            break
    return upper

def main():
    args = get_args()
    tb = datetime.now()
    species, subspecies = parse_taxa(args.taxa, args.no_subspecies)
    
    acc_match, acc_unmatch, tax_match, tax_unmatch = search_fasta(args.input, species, subspecies, args.no_subspecies)
    
    results = [[tax_match,"Matched_Records_Taxon_Names.log"],
               [tax_unmatch,"Unmatched_Records_Taxon_Names.log"],
               [acc_match,"Matched_Records_Accession_Numbers.log"],
               [acc_unmatch,"Unmatched_Records_Accession_Numbers.log"]]
    
    os.chdir(args.out_dir)
    for r in results:
        write_log(r[0],r[1])

    filtered_acc_match = read_accs(results[2][1])
    filtered_acc_unmatch = read_accs(results[3][1])
    
    record_index = index_fasta(args.input)    
    write_fasta(filtered_acc_unmatch,"Unmatched_Taxa",record_index)
    write_fasta(filtered_acc_match,"Matched_Taxa",record_index)
    tf = datetime.now()
    te = tf - tb
    print "\n\nFinished. Total elapsed time: {0} (H:M:S)\n\n".format(te)
    
    
if __name__ == '__main__':
    main()
