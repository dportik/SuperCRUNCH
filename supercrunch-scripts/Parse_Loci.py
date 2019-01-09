'''
SuperCRUNCH: Parse_Loci module

Usage: python Parse_Loci.py -i [fasta file] REQUIRED
			-l [tab-delimited text file with locus information] REQUIRED 
			-t [text file with taxon names to cross-reference] REQUIRED
			-o [path to output directory] REQUIRED
			--no_subspecies (OPTIONAL flag, excludes subspecies names from searches)

    Parse_Loci: The goal of this script is to automate the creation of locus/gene-specific fasta files 
    from a larger fasta file of sequences downloaded from GenBank. The user supplies a 
    GenBank fasta file (ie raw sequence database), a text file containing locus/gene labels,
    abbreviations, and description search terms, and a text file containing a list of taxon 
    names to cross-reference. For a sequence to be added to the locus-specific fasta file,
    it must match at least one of the locus/gene search terms and also have a taxon name 
    that is present in the taxon names text file. Taxon names can contain a mix of species 
    (binomial name) and subspecies (trinomial name) labels. The user supplies a text file 
    containing a list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should not be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record.
    All search terms are converted to uppercase for searching, as well as the description line of 
    each record that is searched, so terms are NOT case-sensitive. The locus-specific 
    fasta files are written to the output directory specified. A log file summarizing the 
    number of records written per locus is also written to the output directory, and is 
    called 'Loci_Record_Counts.log'.
    
[-l] The input locus/gene information file should be tab-delimited with three columns. The
    first column is used to label the output fasta file, and the entry should contain no 
    spaces or special characters. The second column contains gene abbreviations. There can
    be multiple terms as long as they are separated by a semi-colon. The third column 
    contains description search terms, and can contain multiple terms separated by a 
    semi-colon. Spaces can be used in the description search terms, but spaces are not useful 
    to include in abbreviation search terms. All search terms are converted to uppercase for
    searching, so terms are NOT case-sensitive.
    
    Example of file structure for locus information:
    
    [Locus/Gene Label]	[Abbreviation Search Terms]		[Description Search Terms]
    CMOS				CMOS;C-MOS						oocyte maturation factor
	CXCR4				CXCR4							chemokine C-X-C motif receptor 4;C-X-C chemokine receptor type 4
	DLL1				DLL1;DLL						distal-less
	DNAH3				DNAH3							dynein axonemal heavy chain 3;dynein axonemal heavy chain 3s
 
[-t] The input taxon name file should simply contain a list of binomial names, with genus 
 	and species separated by a space, on each line. All taxon names are converted to 
 	uppercase for searching, so names are NOT case-sensitive.
 
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
	
	Loci_Record_Counts.log - Tab-delimited text file summarizing the number of records 
                             written per locus.
							
		Example output contents:
		
		Locus_Name	Records_Written									
		CMOS	9402
		CXCR4	273
		DLL1	257
		DNAH3	1063
	
	[locus name].fasta - The locus-specific fasta file containing all records that match
                         at least one locus/gene search term and contain a valid taxon name
                         (that is present in the taxon information file supplied). A fasta
                         file is created for every locus included in the locus information
                         file supplied.

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
import shutil
import time
import numpy as np
from Bio import SeqIO
from datetime import datetime

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Parse_Loci: The goal of this script is to automate the creation of locus/gene-specific fasta files 
    from a larger fasta file of sequences downloaded from GenBank. The user supplies a 
    GenBank fasta file (ie raw sequence database), a text file containing locus/gene labels,
    abbreviations, and description search terms, and a text file containing a list of taxon 
    names to cross-reference. For a sequence to be added to the locus-specific fasta file,
    it must match at least one of the locus/gene search terms and also have a taxon name 
    that is present in the taxon names text file. Taxon names can contain a mix of species 
    (binomial name) and subspecies (trinomial name) labels. The user supplies a text file 
    containing a list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should not be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record.
    All search terms are converted to uppercase for searching, as well as the description line of 
    each record that is searched, so terms are NOT case-sensitive. The locus-specific 
    fasta files are written to the output directory specified. A log file summarizing the 
    number of records written per locus is also written to the output directory, and is 
    called 'Loci_Record_Counts.log'.
    
    DEPENDENCIES: Python: BioPython, NumPy.
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input", required=True, help="REQUIRED: The full path to a fasta file of GenBank sequence data.")
    parser.add_argument("-l", "--loci", required=True, help="REQUIRED: The full path to a text file containing loci information to search for within the fasta file.")
    parser.add_argument("-t", "--taxa", required=True, help="REQUIRED: The full path to a text file containing all taxon names to cross-reference in the fasta file.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: The full path to an existing directory to write output files.")
    parser.add_argument("--no_subspecies", required=False, action='store_true', help="OPTIONAL: Ignore any subspecies labels in both the name database and record searches (only search binomial names).")
    return parser.parse_args()

def parse_loci_terms(f):
    '''
    Input file is tab delimited, three columns - [Locus ID] [Locus abbreviation] [Locus string]
    Columns two and three can contain multiple search items,
    but if so must be separated by semicolon character (;).
    [Locus ID] - Name or abbreviation containing ONLY alphanumeric characters, used to label output files.
    [Locus abbreviation] - Locus abbreviation, can contain any characters but no spaces.
    [Locus string] - word description of locus, can contain any characters including spaces.

    Example with locus CMOS:
    [Locus ID] - CMOS (but not C-MOS or oocyte maturation factor)
    [Locus abbreviation] - CMOS;C-MOS (note two search terms)
    [Locus string] - oocyte maturation factor; oocyte-maturation factor (note two search terms)
    '''
    print "\n\n\nParsing locus information from {}:".format(f)
    loci_info = []
    with open(f, 'r') as fh_f:
        for line in fh_f:
            #if line not blank, split by tab, replace all quotes and remove whitespace
            cols = [l.replace('\"','').strip().upper() for l in line.split('\t') if line.strip()]
            #for sublist in list, split by semi-colon
            split_cols = [c.split(";") for c in cols]
            #should produce the following from the example:
            #example - CO1	CO1;COX1;COI	"cytochrome oxidase subunit 1;oxidase subunit I;oxidase subunit 1;mitochondrion, complete genome"
            #list created - [ ['CO1'], ['CO1', 'COX1', 'COI'], [CYTOCHROME OXIDASE SUBUNIT 1', 'OXIDASE SUBUNIT I', 'OXIDASE SUBUNIT 1', 'MITOCHONDRION, COMPLETE GENOME'] ]
            loci_info.append(split_cols)
    print "\tFound {} loci to search.\n".format(len(loci_info))
    return loci_info

def parse_taxa(f, no_subspecies):
    '''
    Retrieve taxon names from user supplied taxon names file (f).
    Will load species and subspecies names, and account for 
    species names that are only represented in subspecies labels.
    Summarizes number of each depending on whether subspecies
    are included or excluded (no_subspecies), and regardless will return 
    separate lists for species and subspecies, with all names in uppercase.
    '''
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

def lines_fasta_v1(f):
    '''
    Simply count number of lines in a file (f).
    Turned out to be a faster solution.
    Took 6m 16s for a 16.2 GB file with 1,619,723 lines.
    '''
    print "\nAssessing file size of {}:".format(f)
    tb = datetime.now()
    with open(f) as fh_fasta:
        count = sum(1 for line in fh_fasta)
        print "\tFasta file contains {} lines".format(count)
    tf = datetime.now()
    te = tf - tb
    print "\tTotal time to count lines: {0} (H:M:S)\n".format(te)
    return count

def lines_fasta_v2(f):
    '''
    Simply count number of lines in a file (f).
    Turned out to be a slower solution,
    and is no longer used in the workflow.
    Took 9m 20s for a 16.2 GB file with 1,619,723 lines.
    '''
    print "\nAssessing file size of {}:".format(f)
    tb = datetime.now()
    with open(f, 'r') as fh_fasta:
        lcnt = int(0)
        for line in fh_fasta:
            lcnt += 1
    print "\tFasta file contains {} lines".format(lcnt)
    tf = datetime.now()
    te = tf - tb
    print "\tTotal time to count fasta file (v2): {0} (H:M:S)\n".format(te)
	

def index_fasta(f):
    '''
    Use SeqIO to index the fasta file (f), which may be too big 
    to parse into a list. Returns a dictionary-like structure.
    '''
    tb = datetime.now()
    print "\nIndexing fasta file {}, this could take some time...".format(f)
    record_index = SeqIO.index(f, "fasta")
    tf = datetime.now()
    te = tf - tb
    print "\tTotal time to index fasta file: {0} (H:M:S)\n\n".format(te)

    return record_index
    
def iter_loci(loci_info, f, record_index, out_dir, high, species, subspecies, no_subspecies):
    '''
    Iterate over locus list (loci_info) to search fasta file (f) for 
    matches to taxon names and loci using function 'search_fasta'. Take
    collected information and write fasta file for locus using function
    'write_fasta'.
    '''
    os.chdir(out_dir)
    with open("Loci_Record_Counts.log",'a') as fh_log:
        fh_log.write("Locus_Name\tRecords_Written\n")
    for l in loci_info:
        print "\n\tSearching for {}:".format(l[0][0])
        tb = datetime.now()        
        search_set = search_fasta(f, l[1], l[2], high, species, subspecies, no_subspecies)
        tf = datetime.now()
        te = tf - tb
        print "\t\tSearch time for {0}: {1} (H:M:S)".format(l[0][0],te)
        write_fasta(search_set,l[0][0],record_index)
        with open("Loci_Record_Counts.log",'a') as fh_log:
            fh_log.write("{}\t{}\n".format(l[0][0], len(search_set)) )
            
def get_accession(line):
    '''
    Extracts accession number from a description line in fasta file 
    Example:
    >AF259258.1 Crotalus atrox 12S ribosomal RNA gene, partial sequence; mitochondrial gene for mitochondrial product
    Will yield:
    AF259258.1
    '''
    acc = [l.strip('>') for l in line.split()][0]
    return acc

def taxon_match_sp(taxon_sp,species):
    '''
    Simple function to check if the species (binomial) name
    is present in the species (binomial) list obtained from
    the taxon file. Both will be supplied as uppercase. 
    Returns True or False.
    '''
    match = False
    if taxon_sp in species:
        match = True
    return match

def taxon_match_ssp(taxon_ssp,subspecies):
    '''
    Simple function to check if the subspecies (trinomial) name
    is present in the subspecies (trinomial) list obtained from
    the taxon file. Both will be supplied as uppercase.
    Returns True or False.
    '''
    match = False
    if taxon_ssp in subspecies:
        match = True
    return match

def get_taxon(line):
    '''
    Retrieve the taxon name from '>' line in a fasta file.
    Will fetch the species (binomial) and subspecies (trinomial)
    names and return both. Input line was converted to uppercase,
    so names are also in uppercase.
    '''
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if line.split() >= int(3)][1:3]
    taxon_sp = " ".join(parts1)
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if line.split() >= int(4)][1:4]
    taxon_ssp = " ".join(parts2)
    return taxon_sp, taxon_ssp

def taxon_match_on_line(line, species, subspecies, no_subspecies):
    '''
    Check if taxon name in fasta description line (line) is in taxon lists
    (species, subspecies), depending on whether subspecies are desired 
    or not (no_subspecies = True or False).
    '''
    match = False
    taxon_sp, taxon_ssp = get_taxon(line)
    #Statement for testing if species and subspecies names
    #occur in the line. This should handle all cases correctly.
    #First set of statements are when subspecies are desired:
    if no_subspecies is False:
        if taxon_match_sp(taxon_sp,species) is True:
            if taxon_match_ssp(taxon_ssp,subspecies) is True:
                match = True
                #print "{} found in species and {} found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize()) 
            else:
                match = True
                #print "{} found in species but {} NOT found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize()) 
        elif taxon_match_sp(taxon_sp,species) is False:
            if taxon_match_ssp(taxon_ssp,subspecies) is True:
                match = True
                #print "{} NOT found in species but {} found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize())
            else:
                match = False
                #print "{} NOT found in species and {} NOT found in subspecies".format(taxon_sp.capitalize(),taxon_ssp.capitalize())
    #Second set of statements are when subspecies are NOT desired
    if no_subspecies is True:
        if taxon_match_sp(taxon_sp,species) is True:
            match = True
            #print "{} found in species".format(taxon_sp.capitalize()) 
        elif taxon_match_sp(taxon_sp,species) is False:
            match = False
            #print "{} NOT found in species".format(taxon_sp.capitalize())
    return match

def locus_searchv1(line,terms1):
    '''
    ***Works but is VERY inefficient. No longer incorporated in the main workflow.
    Search for gene abbreviations within parentheses in the
    description line of a fasta file.

    For example in line:
    >AY737396.1 Oligosoma grande cytochrome b (CYTB) gene, partial cds;

    The 'CYTB' term from the parentheses is matched against all abbreviations
    provided in column two of the loci information file. All parentheses are
    detected in this way and searched for.
    '''
    match = False
    if (len(line.split("(")) > int(1)):
    	#First split by ( -> SOMETHING (TARGET) SOMETHING ELSE = [SOMETHING ], [TARGET) SOMETHING ELSE]
    	#Then split by ) -> [SOMETHING ], [ [TARGET], [ SOMETHING ELSE] ]
    	#Use list comprehension to break up all parentheses, all targets will be in list[1:][0]
        ps = [l.split(")") for l in line.split("(")]
        for t in terms1:
            for p in ps[1:]:
                if t == p[0]:
                    match = True
    return match

def locus_searchv2(line,terms1,terms2):
    '''
    Search for gene descriptions (terms1; list) and abbreviations (terms2; list) 
    in the description line (line) of a fasta file. Works much faster and
    also finds all records that the above parentheses search function will find.
    '''
    match = False
    #split line by whitespace and remove following characters -> , ; ( )
    p = [l.replace(",",'').replace(";",'').replace(":",'').replace("(",'').replace(")",'') for l in line.split()]
    for t in terms1:
        if t in p:
            match = True
    for t in terms2:
        if t in line:
            match = True
    return match

def find_interval(count):
    '''
    Find where a number fits in a distribution. Will find 
    when the number is exceeded by a value in the distribution
    and returns that value. Needed for tracking progress on
    big lists (file lines, accession numbers, etc).
    '''
    vala = np.array([100,1000,10000])
    valb = np.arange(100000,1000000000000,100000)
    vals = np.append(vala,valb)
    for i in vals:
        if count <= i:
            high = i
            break
    return high

def search_fasta(f, terms1, terms2, high, species, subspecies, no_subspecies):
    '''
    Iterate over lines in fasta file (f), finds record lines (with >),
    checks if locus information (terms1, terms2; both lists) are found in description 
    line. If so, then checks if taxon name is in taxon names lists (species,
    subspecies), depending on subspecies option (no_subspecies = True or False). 
    If everything matches, will add accession number to the set. Returns set 
    of accession numbers, which is used for writing a new fasta file for the locus in 
    function 'write_fasta' below.
    '''
    search_set = set()
    with open(f, 'r') as fh_fasta:
        print "\t\tSearch progress:"
        lcnt = int(0)
        for line in fh_fasta:
            #setup basic line counter to show progress on large files
            lcnt += 1
            interval = high // 20
            if lcnt % interval == 0:
                print "\t\t\tLine {}...".format(lcnt)
            #find records by >
            if line.startswith('>'):
                if '|' in line:
                    line.replace("|"," ")
                #convert line to uppercase and remove line break
                line = line.upper().strip()
                #search locus terms with function, will return T or F
                if locus_searchv2(line,terms1,terms2) is True:
                    #check if taxon is in databse with function, will return T or F
                    if taxon_match_on_line(line, species, subspecies, no_subspecies) is True:
                        #if above is all passed, add accession to set
                        #extract accession number with function
                        acc = get_accession(line)
                        #add to set
                        search_set.add(acc)
    return search_set

def write_fasta(search_set,prefix,record_index):
    '''
    Write fasta file using BioPython based on set of accession numbers 
    provided (search_set) and dictionary from indexed fasta file (record_index).
    '''
    out_fasta = "{}.fasta".format(prefix)
    acc_list = sorted(search_set)
    with open(out_fasta, 'a') as fh_out:
        for acc in acc_list:
            fh_out.write(">{0}\n{1}\n".format(record_index[acc].description,record_index[acc].seq))
    print "\t\tWrote {0} records to {1}\n".format(len(search_set),out_fasta)

def main():
    '''
    Get arguments
    Parse loci information file and return list
    Parse taxon names file and return list
    Obtain line count of input fasta file
    Find suitable number for line progress report
    Index the input fasta file (can be a long step)
    Iterate over loci in loci list to perform searches and write output fasta files
    '''
    args = get_args()
    tb = datetime.now()
    loci_info = parse_loci_terms(args.loci)
    species, subspecies = parse_taxa(args.taxa, args.no_subspecies)
    count = lines_fasta_v1(args.input)
    high = find_interval(count)
    record_index = index_fasta(args.input)
    iter_loci(loci_info, args.input, record_index, args.out_dir, high, species, subspecies, args.no_subspecies)
    tf = datetime.now()
    te = tf - tb
    print "Total time to parse all loci: {0} (H:M:S)\n\n".format(te)
    
    
if __name__ == '__main__':
    main()
