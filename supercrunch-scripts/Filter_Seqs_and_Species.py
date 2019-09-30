'''
SuperCRUNCH: Filter_Seqs_and_Species module
                             
    Filter_Seqs_and_Species: A tool for filtering sequences from a fasta file to 
    select a single representative sequence per taxon or select all sequences per 
    taxon. Taxon names in the decription line of sequence records are matched
    to those from a supplied taxon names file in a similar way to Parse_Loci.py
    (through an SQL database constructed from the sequence descriptions). 
    Taxon names can contain a mix of species (two-part) and subspecies 
    (three-part) names. The --no_subspecies flag can be used to only 
    include binomial names in searches.

    The main decision for filtering is the -s flag, with options 'oneseq' or 'allseqs'.
    If 'oneseq' is selected, a single representative sequence will be selected for
    each taxon. This will lead to a supermatrix-style data set. If 'allseqs' is 
    selected, all sequences will be selected for each taxon. This will allow
    population-level (intraspecific) sampling. In both cases, the sequences must 
    be longer than the base pair length requirement supplied using the -m flag. 
    If no sequences are long enough, NO sequence is selected for that taxon. 
    
    The remaining filters generally apply to the 'oneseq' option. When sequences
    are found for a taxon, they are first sorted by length (longest first). 
    A 'best' sequence is selected based on the -f flag. The 'length' option takes the 
    longest sequence. The 'translate' option takes the longest sequence
    that passes translation (no internal stop codons). For this option the 
    translation table should be specified with the --table flag. All NCBI 
    translation table options are available. If the --table flag is omitted, 
    the default will be to use the Standard translation table.    
    The --randomize flag can be provided to shuffle the sequences in random
    order, rather than by length. If this is used with -f length, sequence 
    selection is essentially randomized (but minimum length filter must be passed). 

    If the --vouchered flag is provided, only sequences which have a 
    voucher tag (Voucher_[ID]) will be retained. This can be useful for creating 
    a population-level data set for which particular specimens are known to
    have been sequenced for several loci. It will allow them to be concatenated 
    into a population-level supermatrix downstream. 

    Finally, if your dataset contains a mix of loci that need to be filtered 
    differently (e.g., coding mtDNA vs. coding nucDNA vs. noncoding), the 
    --onlyinclude flag can be used. This requires the full path to a text file 
    containing the names of the fasta files (one file name per line) to be 
    processed for a particular run. 

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.

-------------------------
Compatible with Python 2.7 & 3.7
Python packages required:
	-BioPython
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
import shutil
import itertools
import operator
import random
import sqlite3
from datetime import datetime
from random import shuffle
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Filter_Seqs_and_Species: A tool for filtering sequences from a fasta file to 
    select a single representative sequence per taxon or select all sequences per 
    taxon. Taxon names in the decription line of sequence records are matched
    to those from a supplied taxon names file in a similar way to Parse_Loci.py
    (through an SQL database constructed from the sequence descriptions). 
    Taxon names can contain a mix of species (two-part) and subspecies 
    (three-part) names. The --no_subspecies flag can be used to only 
    include binomial names in searches.

    The main decision for filtering is the -s flag, with options 'oneseq' or 'allseqs'.
    If 'oneseq' is selected, a single representative sequence will be selected for
    each taxon. This will lead to a supermatrix-style data set. If 'allseqs' is 
    selected, all sequences will be selected for each taxon. This will allow
    population-level (intraspecific) sampling. In both cases, the sequences must 
    be longer than the base pair length requirement supplied using the -m flag. 
    If no sequences are long enough, NO sequence is selected for that taxon. 
    
    The remaining filters generally apply to the 'oneseq' option. When sequences
    are found for a taxon, they are first sorted by length (longest first). 
    A 'best' sequence is selected based on the -f flag. The 'length' option takes the 
    longest sequence. The 'translate' option takes the longest sequence
    that passes translation (no internal stop codons). For this option the 
    translation table should be specified with the --table flag. All NCBI 
    translation table options are available. If the --table flag is omitted, 
    the default will be to use the Standard translation table.    
    The --randomize flag can be provided to shuffle the sequences in random
    order, rather than by length. If this is used with -f length, sequence 
    selection is essentially randomized (but minimum length filter must be passed). 

    If the --vouchered flag is provided, only sequences which have a 
    voucher tag (Voucher_[ID]) will be retained. This can be useful for creating 
    a population-level data set for which particular specimens are known to
    have been sequenced for several loci. It will allow them to be concatenated 
    into a population-level supermatrix downstream. 

    Finally, if your dataset contains a mix of loci that need to be filtered 
    differently (e.g., coding mtDNA vs. coding nucDNA vs. noncoding), the 
    --onlyinclude flag can be used. This requires the full path to a text file 
    containing the names of the fasta files (one file name per line) to be 
    processed for a particular run. 

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa'. The 
    NAME portion should not contain any periods or spaces, but can contain 
    underscores. Output files are labeled using a prefix identical to NAME.
    		
    DEPENDENCIES: Python: BioPython.
	-----------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the fasta files to filter by species and sequence.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-s", "--seq_selection",
                            required=True,
                            choices=["oneseq", "allseqs"],
                            help="REQUIRED: Select one representative sequence per taxon "
                            "or select all sequences available per taxon.")
    
    parser.add_argument("-f", "--seq_filter",
                            required=True,
                            choices=["translate", "length"],
                            help="REQUIRED: Strategy for filtering sequence data, "
                            "particularly for selecting one representative sequence "
                            "using the -s 'oneseq' option.")
    
    parser.add_argument("-m", "--min_length",
                            required=True,
                            help="REQUIRED: Integer of the minimum number of base pairs "
                            "required to keep a sequence (ex. 150).")
    
    parser.add_argument("-t", "--taxa",
                            required=True,
                            help="REQUIRED: The full path to a text file containing all "
                            "taxon names to cross-reference in the fasta file(s).")
    
    parser.add_argument("--no_subspecies",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Ignore any subspecies labels in both the name "
                            "database and record searches. Only searches binomial "
                            "(two-part) names).")
    
    parser.add_argument("--randomize",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: For taxa with multiple sequences, shuffle order "
                            "randomly. Overrides sorting by length for sequence selection step.")
    
    parser.add_argument("--vouchered",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Select only sequences that contain a Voucher_[ID] "
                            "tag in the description line (inserted during Parse_Loci.py).")
    
    parser.add_argument("--table",
                            default="standard",
                            choices=["standard","vertmtdna","invertmtdna","yeastmtdna","plastid",
                                         "1","2","3","4","5","6","7","8","9","10","11","12","13",
                                         "14","15","16","17","18","19","20","21","22","23","24",
                                         "25","26","27","28","29","30","31"],
                            help="REQUIRED for -f translate: Specifies translation table.")
    
    parser.add_argument("--onlyinclude",
                            required=False,
                            help="OPTIONAL: The full path to a text file containing the "
                            "names of the fasta files to process for this run.")
    
    parser.add_argument("--quiet",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: Show less output while running (useful when filtering "
                            "many, many loci).")
    
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

def get_taxon(line):
    """
    Retrieve the taxon name from '>' line in a fasta file.
    Will fetch the species (binomial) and subspecies (trinomial)
    names and return both. Input line gets converted to uppercase,
    so names are also in uppercase.
    """
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if line.split() >= int(3)][1:3]
    taxon_sp = " ".join(parts1)
    
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if line.split() >= int(4)][1:4]
    taxon_ssp = " ".join(parts2)
    
    return taxon_sp, taxon_ssp


def translation_val(tcode):
    """
    Convert shortcut translation table terms into the true NCBI
    nomenclature using a dictionary structure. Term needed for 
    translation with biopython.
    """
    tdict = {"standard":"Standard", "vertmtdna":"Vertebrate Mitochondrial",
                 "invertmtdna":"Invertebrate Mitochondrial", "yeastmtdna":"Yeast Mitochondrial",
                 "plastid":"11", "1":1, "2":2, "3":3, "4":4, "5":5, "6":6, "7":7, "8":8, "9":9,
                 "10":10, "11":11, "12":12, "13":13, "14":14, "15":15, "16":16, "17":17, "18":18,
                 "19":19, "20":20, "21":21, "22":22, "23":23, "24":24, "25":25, "26":26, "27":27,
                 "28":28, "29":29, "30":30, "31":31}
    tval = tdict[tcode]
    
    return tval

def check_translation(newseq, tcode):
    """
    Function to check translation of a sequence (newseq)
    based on a specified translation table (tcode). Ultimately
    returns True or False if sequence was successfully translated.
    """
    #adjust starting sequence length to be divisible by three for translation check
    n = Seq("N")
    if len(newseq) % 3 == 0:
        pass
    elif len(newseq) %3 != 0:
        newseq = newseq + n
        if len(newseq) % 3 == 0:
            pass
        elif len(newseq) % 3 != 0:
            newseq = newseq + n

    #Translate based on table selected by user
    #translate to all forward frames, constantly adjusting
    #to keep multiples of 3
    frame1 = newseq[:].translate(table=tcode)
    f2 = newseq+n
    frame2 = f2[1:].translate(table=tcode)
    f3 = f2+n
    frame3 = f3[2:].translate(table=tcode)
    #reverse complement and translate to all frames, constantly adjusting
    #to keep multiples of 3
    rcseq = newseq.reverse_complement()
    rcframe1 = rcseq[:].translate(table=tcode)
    rf2 = rcseq+n
    rcframe2 = rf2[1:].translate(table=tcode)
    rf3 = rf2+n
    rcframe3 = rf3[2:].translate(table=tcode)

    #count stop codons present in each frame
    f1stops = frame1.count('*')
    f2stops = frame2.count('*')
    f3stops = frame3.count('*')
    rcf1stops = rcframe1.count('*')
    rcf2stops = rcframe1.count('*')
    rcf3stops = rcframe1.count('*')

    #create sublists with relevant information to iterate over
    #for each frame - [stop codon count, sequence object in that frame, frame string label]
    frames_list = [[f1stops, frame1, "Frame 1"],
                  [f2stops, frame2, "Frame 2"],
                  [f3stops, frame3, "Frame 3"],
                  [rcf1stops, rcframe1, "RC Frame 1"],
                  [rcf2stops, rcframe2, "RC Frame 2"],
                  [rcf3stops, rcframe3, "RC Frame 3"]]
        
    #create multi-conditional statement
    #first, we search all frames to see if they have 0 stop codons,
    #if so we use that sequence frame
    #second, we search all frames to see if a single stop codon occurs in
    #the final two positions, if so we count as pass
    translated = False
    #sublist will have following contents [stop codon count, sequence in frame, frame label]
    for sublist in frames_list:
        if translated is False:        
            #check if frame has no stop codons
            if sublist[0] == int(0):
                translated = True

    if translated is False:
        for sublist in frames_list:
            if translated is False:
                #check if frame has one stop codon and it is in the last three
                #positions of the sequence
                if sublist[0] == int(1) and '*' in sublist[1][-4:]:
                    translated = True
        
    return translated

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

def parse_fasta_record(line, species, subspecies, no_subspecies):
    """
    Deconstruct a record line from a fasta file into
    several elements which will be added to the SQL
    database. Returns an accession number, species label,
    subspecies label, and T/F for whether voucher tag is
    present. Name labels are constructed based on several
    arguments provided.
    """
    #get accession number from line
    accession = [l.strip('>') for l in line.split()][0]

    #look for voucher flag inserted by Parse_Loci.py:
    if any(i.startswith('Voucher_') for i in line.split()):
        voucher = 'yes'
    else:
        voucher = 'no'
    
    #convert line to uppercase for other strings
    line = line.upper()
    
    #get the 'species' name - first two strings following
    #the accession number
    sp_parts = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                  if len(line.split()) >= int(3)][1:3]
    #construct name and test if it is in the taxon list (species),
    #generates actual name or "NA"
    sp_name = test_name(species, sp_parts, True)

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
        
    return accession, sp_name, ssp_name, voucher

def filter_seqs(locus, cur, taxa, seq_selection, seq_filter,
                    length, randomize, no_subspecies, vouchered, quiet):
    """
    Function to select sequence(s). Heavily annotated
    to explain decision tree, options, and variables
    returned.
    """
    if quiet:
        print("\tSelecting and writing sequences.")
        
    #create a list that will simply contain all the accession numbers for
    #every sequence and every taxon for this locus (e.g. this fasta file)
    all_acc_list = []

    #create a list that will contain ONLY the accession numbers for
    #the sequences that passed filters and will be written to the new
    #filtered fasta file
    writing_acc_list = []

    #create a list with more complete information for the sequences
    #that passed filters and will be written to the new filtered fasta
    #file. each sequence gets a sublist with contents:
    #[taxon name, accessoin, length, vouchered?, translatable?, number of other seqs available]
    all_seq_info = []

    #create a list that contains information about the sequences available for
    #each taxon. each taxon gets a sublist with contents:
    #[taxon name, accession string]
    #the accession string is all the accession numbers for this taxon
    #joined by ', ' - which results in 'Acc1, Acc2, Acc3, ...'
    taxon_accs = []
    
    for t in taxa:        
        #If taxon is a binomial name, exclude the sspname from the search.
        #In creating the records, if the no_subspecies flag was included then
        #ALL the records got an "NA" for the sspname column. If the no_subspecies
        #flag was included then only records with valid subspecies name got a
        #value in the sspname column, and those without valid subspecies names
        #got an "NA" assigned. Therefore, if a taxon name is binomial, we only
        #want records that have that exatc name in spname and an "NA" in sspname.
        #Relax brain, this actually works.
        if len(t.split()) == 2:
            #get all the records that match this taxon name for this locus
            cur.execute("SELECT * FROM records WHERE spname = '{0}' AND sspname = 'NA'"
                            .format(t))
            
            allresults = cur.fetchall()
            
            #get all sequences that match for this name type and this locus and with minimum length required
            cur.execute("SELECT * FROM records WHERE spname = '{0}' AND sspname = 'NA' AND length >= {1}"
                            .format(t, length))
            results = cur.fetchall()

        #If the taxon name is a 'subspecies' name, then we only want records
        #that match that name in the sspname column
        elif len(t.split()) == 3:
            #get all the records that match this taxon name for this locus
            cur.execute("SELECT * FROM records WHERE (spname = '{0}' OR sspname = '{0}')"
                            .format(t))
            allresults = cur.fetchall()
            
            #get all sequences that match for this name type and this locus and with minimum length required
            cur.execute("SELECT * FROM records WHERE (spname = '{1}' OR sspname = '{1}') AND length >= {2}"
                            .format(locus, t, length))
            results = cur.fetchall()
            
            
        #Note that if we do not take the above approach we can end up with
        #duplicate sequence records in our fasta file, because all subspecies
        #will also match to the species name and will be included under that
        #species name if we are using the allseqs option.
        
        #make a sorted list of unique accession numbers
        availableseqs = sorted(set([r['accession'] for r in allresults]))

        #add all accession numbers to the unfiltered accession list
        #which can be used to download all seqs for this fasta file
        [all_acc_list.append(r['accession']) for r in allresults]

        if availableseqs:
            #create a string of the accs for this taxon for one of the
            #output files which looks like: Genus species Acc1, Acc2, Acc3, ...
            acc_join = ", ".join(availableseqs)
            taxon_accs.append([t.capitalize(), acc_join])

        
        resultlist = [[r['accession'], r['length'], r['voucher'], r['translation'],
                           r['spname'].capitalize(), r['sspname'].capitalize()] for r in results]

        #make sure we have results
        if resultlist:
            if not quiet:
                print("\n\n\t{}:".format(t.capitalize()))
                print("\n\t\tNumber of unfiltered sequences available: {}".format(len(availableseqs)))

            #if only vouchered seqs desired, use list comprehension
            #to find entries that have yes in voucher column
            if vouchered:
                resultlist = [x for x in resultlist if x[2] == 'yes']

            #check again because filtering by vouchered can produce nothing
            if resultlist:                    
                #sort by sequence length, longest seqs first
                resultlist.sort(key=operator.itemgetter(1), reverse=True)

                #sort by longest sequence or randomly depending on arg randomize
                if randomize:
                    random.shuffle(resultlist)

                #now that the sequence set is pre-filtered, next we want to select
                #a sequence or sequences based on two remaining arguments:
                #seq_filter and seq_selection
                #begin with seq_filter and then move to seq_selection within
                
                #option for choosing the longest sequence available
                if seq_filter == "length":
                    if not quiet:
                        print("\t\tNumber of filtered sequences available: {}".format(len(resultlist)))
                        print("\n\t\tFiltered sequences: [Accession, Sequence Length]")
                        for r in resultlist:
                            print("\t\t\t\t[{}, {}]".format(r[0], r[1]))
                    
                    #option for taking only the longest sequence
                    if seq_selection == "oneseq":
                        writing_acc_list.append(resultlist[0][0])
                        if not quiet:
                            print("\n\t\tSelecting sequence: {}".format(resultlist[0][0]))
                        #create list that contains taxon name, accession, length, vouchered, translation, availableseqs
                        all_seq_info.append([t.capitalize(), resultlist[0][0],
                                                 resultlist[0][1], resultlist[0][2],
                                                 resultlist[0][3], len(availableseqs)])

                    #option for taking all the available sequences
                    elif seq_selection == "allseqs":
                        for x in resultlist:
                            writing_acc_list.append(x[0])
                            all_seq_info.append([t.capitalize(), x[0],
                                                     x[1], x[2],
                                                     x[3], len(availableseqs)])
                        if not quiet:
                            print("\n\t\tSelecting all sequences above.")

                #option for choosing the longest translatable sequence available
                elif seq_filter == "translate":
                    #filter the resultlist to obtain sequences that pass translation
                    transresultlist = [x for x in resultlist if x[3] == 'yes']
                    if not quiet:
                        print("\t\tNumber of sequences passing translation: {}".format(len(transresultlist)))
                    #the above action can produce an empty list (e.g., no sequences passed
                    #translation), so before moving on to seq_selection make sure there
                    #are actually sequences here
                    if transresultlist:
                        if not quiet:
                            print("\n\t\tFiltered sequences: [Accession, Sequence Length]")
                            for r in transresultlist:
                                print("\t\t\t\t[{}, {}]".format(r[0], r[1]))
                        
                        #option for taking only the longest translatable sequence
                        if seq_selection == "oneseq":
                            writing_acc_list.append(transresultlist[0][0])
                            if not quiet:
                                print("\n\t\tSelecting sequence: {}".format(transresultlist[0][0]))
                            #add sublist of [taxon name, accession, length, vouchered, translation, availableseqs] to larger list
                            all_seq_info.append([t.capitalize(), transresultlist[0][0],
                                                     transresultlist[0][1], transresultlist[0][2],
                                                     transresultlist[0][3], len(availableseqs)])

                        #option for taking all translatable sequences
                        elif seq_selection == "allseqs":
                            accs = sorted([r[0][0] for r in transresultlist])
                            for x in transresultlist:
                                writing_acc_list.append(x[0])
                                #add sublist of [taxon name, accession, length, vouchered, translation, availableseqs] to larger list
                                all_seq_info.append([t.capitalize(), x[0],
                                                         x[1], x[2],
                                                         x[3], len(availableseqs)])
                            if not quiet:
                                print("\n\t\tSelecting all sequences above.")
                            
                    #if no sequences passed translation, handle here 
                    else:
                        if not quiet:
                            print("\n\t\t***All sequences failed translation!")
                        
                        #if the goal is to select one sequence per taxon, we still
                        #want a sequence for this taxon even if none of them passed
                        #translation. in this case, we will just take the longest sequence
                        if seq_selection == "oneseq":
                            writing_acc_list.append(resultlist[0][0])
                            if not quiet:
                                print("\t\tSelecting sequence: {}".format(resultlist[0][0]))
                            #create list that contains taxon name, accession, length, vouchered, translation, availableseqs
                            all_seq_info.append([t.capitalize(), resultlist[0][0],
                                                     resultlist[0][1], resultlist[0][2],
                                                     resultlist[0][3], len(availableseqs)])
                                                        
            #this is why argument vouchered can be a bad idea!
            #look at how it overrides the main decision tree to
            #easily eliminate all available sequences!
            else:
                if not quiet:
                    print("\n\t\tNo vouchered seqs were found!")
    
    return all_acc_list, writing_acc_list, all_seq_info, taxon_accs
    
def write_fasta(f, locus, seq_selection, writing_acc_list, fdir):
    """
    Function to write all sequences to a new fasta
    file for a locus, based on the accession numbers
    provided. Uses biopython to index fasta file and 
    exploit dictionary-like structure to quickly write
    seqs using the accession numbers. Moves output to the
    appropriate output directory.
    """
    recs = SeqIO.index(f, "fasta")
    outname = "{}_{}.fasta".format(locus, seq_selection)
    
    with open(outname, 'a') as fh:
        for a in writing_acc_list:
            fh.write((recs[a]).format("fasta"))
            
    if os.stat(outname).st_size == 0:
        print("\nNo records written for {}, removing empty file.\n".format(outname))
        os.remove(outname)
        
    else:
        shutil.move(outname, fdir)
        
def write_species_log(locus, all_seq_info, ldir):
    """
    Function to write detailed information for all
    sequences written to the filtered fasta file,
    using the all_seq_info list. Moves output to the
    appropriate output directory.
    """
    outname = "{}_species_log.txt".format(locus)
    
    with open(outname, 'a') as fh:
        fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
                     .format("Taxon", "Acession", "SeqLength", "Vouchered", "PassedTranslation", "SeqsAvailable"))
        
        for x in all_seq_info:
            fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
                         .format(x[0], x[1], x[2], x[3], x[4], x[5]))
            
    shutil.move(outname, ldir)
            
def write_species_accs(locus, taxon_accs, adir):
    """
    Function to write a two-column output file which
    contains all taxa and the associated accession numbers.
    Moves output to the appropriate output directory.
    """
    outname = "{}_accession_list_by_species.txt".format(locus)
    with open(outname, 'a') as fh:
        for x in taxon_accs:
            fh.write("{0}\t{1}\n".format(x[0], x[1]))
            
    shutil.move(outname, adir)

def write_batch_entrez_accs(locus, all_acc_list, adir):
    """
    Function to write simple list of all accession
    numbers that were in the starting fasta file, 
    for the purpose of being able to reconstruct this
    file by downloading accs using batch entrez. 
    Moves output to the appropriate output directory.
    """
    outname = "{}_accessions_for_BatchEntrez.txt".format(locus)
    with open(outname, 'a') as fh:
        for x in all_acc_list:
            fh.write("{0}\n".format(x))
            
    shutil.move(outname, adir)
    
def make_locus_db(f, species, subspecies, no_subspecies, seq_filter, tcode, quiet):
    """
    Add the entire contents of this fasta file (f)
    to the locus SQL database.
    """
    #create name for this db based on locus name
    locusdb =  "{}.sql.db".format(f.split('/')[-1].split('.')[0])
    curpath = os.getcwd()
    db = os.path.join(curpath, locusdb)
    #double check it doesn't already exist - if it does eliminate
    if os.path.exists(db):
        os.remove(db)
    #establich connection to db
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    #create table in db
    cur.execute("DROP TABLE IF EXISTS records;")
    cur.execute("""CREATE TABLE records (
        accession text NOT NULL,
        spname text NOT NULL,
        sspname text NOT NULL,
        voucher text NOT NULL,
        length int NOT NULL,
        translation text NOT NULL)""")
    sql_add_records = "INSERT INTO records VALUES (?, ?, ?, ?, ?, ?)"

    if quiet:
        print("\tBuilding locus SQL database: {}".format(locusdb))
    elif not quiet:
        print("\nBuilding locus SQL database: {}".format(locusdb))
        
    #load fasta as seq dict using biopython
    locus_recs = SeqIO.index(f, "fasta")
    #iterate over fasta contents and add info to sql db
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                #get main contents using parse_fasta_record() function
                acc, sp, ssp, vouch = parse_fasta_record(line, species, subspecies, no_subspecies)
                #get access to sequence to potentially translate depending on arg supplied
                tempseq = locus_recs[acc].seq
                if seq_filter == "translate":
                    if check_translation(tempseq, tcode):
                        tran = "yes"
                    else:
                        tran = "no"
                else:
                    tran = "NA"
                #add information to sql db
                cur.execute(sql_add_records, (acc, sp, ssp, vouch, len(tempseq), tran))
                
    #save changes to db, close connection, return name of db
    conn.commit()
    conn.close()
    
    return db

def filter_runner(flist, species, subspecies, seq_selection, seq_filter, length,
                      randomize, no_subspecies, vouchered, tcode, quiet,
                      fdir, ldir, adir):
    """
    Function to run main tasks. Connect to the 
    constructed SQL database, get taxon names from
    the empirical fasta files (using the db), then 
    iterate over every fasta file (locus) to select
    the relevant sequences. Writes output files and
    cleans up for each locus.
    """
    print("\n--------------------------------------------------------------------------------------\n")
    print("Beginning sequence filtering for all loci.\n\n\n")
    
    #iterate over loci
    for f in flist:
        b = datetime.now()
        #get locus name from full path to fasta file
        locus = f.split('/')[-1].split('.')[0]
        if not quiet:
            print("\n--------------------------------------------------------------------------------------\n")
        print("Filtering sequences in {}.".format(locus))

        #make an SQL db for this locus
        db = make_locus_db(f, species, subspecies, no_subspecies, seq_filter, tcode, quiet)
        conn = sqlite3.connect(db)
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()

        #get a list of the species and subspecies from this particular file
        #using the sql db
        cur.execute("SELECT * FROM records")
        results = cur.fetchall()
        emp_species = sorted(set([r['spname'] for r in results if r['spname'] != "NA"]))
        emp_subspecies = sorted(set([r['sspname'] for r in results if r['sspname'] != "NA"]))

        #create taxon list to sort through depending on
        #the no_subspecies option
        if no_subspecies:
            taxa = emp_species
        else:
            taxa = sorted(emp_species + emp_subspecies)

        #filter the sequences using the filter_seqs() function
        all_acc_list, writing_acc_list, all_seq_info, taxon_accs = filter_seqs(locus, cur, taxa,
                                                                                   seq_selection, seq_filter,
                                                                                   length, randomize, no_subspecies,
                                                                                   vouchered, quiet)
        #write the filtered fasta file and all log files
        write_fasta(f, locus, seq_selection, writing_acc_list, fdir)
        write_species_log(locus, all_seq_info, ldir)
        write_species_accs(locus, taxon_accs, adir)
        write_batch_entrez_accs(locus, all_acc_list, adir)
        #remove the sql db for this locus
        os.remove(db)
        #show elapsed time on screen
        f = datetime.now()
        e = f - b
        if quiet:
            print("\tFinished. Elapsed time: {0} (H:M:S)\n".format(e))
        elif not quiet:
            print("\n\nFinished. Elapsed time: {0} (H:M:S)\n".format(e))

def make_dirs():
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    curpath = os.getcwd()
    
    rdir = os.path.join(curpath, "Results")
    if not os.path.exists(rdir):
        os.mkdir(rdir)
        
    fdir = os.path.join(curpath, "Results", "Filtered-Fasta-Files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)
        
    ldir = os.path.join(curpath, "Results", "Taxon-Log-Files")
    if not os.path.exists(ldir):
        os.mkdir(ldir)
        
    adir = os.path.join(curpath, "Results", "Accession-Number-Files")
    if not os.path.exists(adir):
        os.mkdir(adir)
        
    return rdir, fdir, ldir, adir

#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    tb = datetime.now()
    
    species, subspecies = parse_taxa(args.taxa)
    
    tcode = translation_val(args.table)
    
    os.chdir(args.indir)
    if not args.onlyinclude:
        flist = sorted([os.path.abspath(f) for f in os.listdir('.') if f.endswith(('.fasta', '.fa'))])
        print("\nFound {:,} fasta files to filter.".format(len(flist)))
        
    elif args.onlyinclude:
        with open(args.onlyinclude, 'r') as fh:
            incl = [l.strip() for l in fh]
        flist = sorted([os.path.abspath(f) for f in os.listdir('.') if f.endswith(('.fasta', '.fa')) and
                            f in incl])
        print("\nFound {:,} fasta files to filter, based on list of files to include.".format(len(flist)))

    os.chdir(args.outdir)
    rdir, fdir, ldir, adir = make_dirs()
    
    filter_runner(flist, species, subspecies, args.seq_selection, args.seq_filter, args.min_length,
                      args.randomize, args.no_subspecies, args.vouchered, tcode, args.quiet,
                      fdir, ldir, adir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nTotal time to filter all files: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
            
if __name__ == '__main__':
    main()
