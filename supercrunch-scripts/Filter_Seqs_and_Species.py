'''
SuperCRUNCH: Filter_Seqs_and_Species module

Usage: python Filter_Seqs_and_Species.py -i [directory containing both files] REQUIRED
                             -f [sequence filtering method (translate,length)] REQUIRED
                             -l [minimum number of base pairs required to pass] REQUIRED
                             -t [text file with taxon names to cross-reference] REQUIRED
                             --no_subspecies (OPTIONAL; exclude subspecies names from searches)
                             --randomize (OPTIONAL; sort seqs randomly rather than by length)
                             --allseqs (OPTIONAL; includes all seqs for all taxa)
                             --table (OPTIONAL; DEFAULT = standard)
                             
    Filter_Seqs_and_Species: A tool for filtering sequences from a fasta file to filter and select a single
    representative sequence per taxon or to filter and select all sequences per taxon. 
    The script examines taxon names in the decription line of sequence records to see if
    they match those from a taxonomic database. Taxon names can contain a mix of species (binomial name) 
    and subspecies (trinomial name) labels. The user supplies a text file containing a
    list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should not be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record. For all
    valid taxa names found, the taxon names are compiled and for each taxon all records are found.

    The sequences available for a given taxon are sorted by length (longest first). A sequence is 
    selected based on the -f flag (options: translate,length), where in the case of 
    coding loci the first sequence in the list that passes a  translation test is selected. 
    If no sequences pass translation the first (longest) sequence is selected. 

    The translation table should be specified with the --table flag. All NCBI translation table options 
    are available, and can be selected using integers or the shortcut terms provided. If the --table 
    flag is omitted, the default will be to use the Standard translation table. 
    
    For the length strategy (-f length) the first (longest) sequence is selected. 
    In addition, all sequences must be longer than the base pair length requirement supplied 
    (-l). If a sequence is not long enough, the search continues down the list of available 
    sequences. If no sequences are long enough, NO sequence is selected for that taxon. 

    Both strategies (translate, length) provide an objective and repeatable method to reconstruct 
    an identical supermatrix every run. Finally, if the --randomize flag is provided, the list of 
    available sequences for a given taxon is shuffled randomly, rather than sorted by length. 
    If this is used in conjunction with -f standard or vmtdna, a translation test and length test 
    are still performed to select a sequence. If this is used with -f length, then the sequence 
    selection is truly randomized (but the length filter must be passed). If there are 
    multiple sequences available for many taxa and loci, the -f length plus --randomize 
    option will provide a nearly infinite number of permutations of the final supermatrix 
    constructed.

    If all sequences passing the filters are desired instead of one representative sequence per
    taxon, the --allseqs flag can be used. This is useful for filtering and retaining intraspecific 
    sequence data, but it will not allow loci to be concatenated into a supermatrix downstream.

    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

    Major Steps:

    1. Finds all species names in fasta file.
    2. Iterates over species names and find accession numbers of all sequence records matching taxon.
    3. Tests if sequence is longer than the user-supplied length requirement, adds to list for taxon.
    4. Sorts sequences by length, then iterates over sequences. If optional flag --randomize is provided,
       the sequences are shuffled into random order before the subsequent steps. This 'randomize'
       feature can be used to generate a nearly infinite number of permutations of the supermatrix
       if there are many sequences available for a given taxon. 
    5a. If standard or vmtdna, sequences undergo a translation test to see if a correct reading
        frame can be identified. If so, it is taken as the sequence for the taxon. If not, moves
        to the next sequence and tests. If no sequences pass translation the first sequence
        (which is the longest) will be used for the taxon.
    5b. If noncoding, the longest sequence is taken for the taxon.


Output files:

    [fasta name]_single_taxon.fasta - The output fasta file which contains a single filtered sequence
                                    per taxon.
                OR
    [fasta name]_all_seqs.fasta - The output fasta file which contains all filtered sequences available
                                    per taxon.


    [fasta name]_species_log.txt - A summary file containing the following columns of information
                                    for each taxon entry:
                                    Taxon  Accession  SeqLength   PassedTranslation  SeqsAvailable

    [fasta name]_accession_list_by_species.txt - A tab-delimited file in which each line starts
                                    with a taxon name and is followed by all accession numbers of
                                    sequences passing the length filter from the fasta file.

    [fasta name]_accession_list_for_Batch_Entrez.txt - a batch entrez style file which simply
                                contains all accession numbers of sequences passing length
                                filter from the fasta file.

-------------------------
For Python 2.7
Python modules required:
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
from datetime import datetime
from random import shuffle
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
    Filter_Seqs_and_Species: A tool for filtering sequences from a fasta file to filter and select a single
    representative sequence per taxon or to filter and select all sequences per taxon. 
    The script examines taxon names in the decription line of sequence records to see if
    they match those from a taxonomic database. Taxon names can contain a mix of species (binomial name) 
    and subspecies (trinomial name) labels. The user supplies a text file containing a
    list of taxon names to cross-reference. All taxon names are converted to uppercase for
    searching, as well as the description line of each record that is searched, so taxon
    names are NOT case-sensitive. If records with valid subspecies names should NOT be included, 
    use the --no_subspecies flag. If this flag is used, only the genus and species will be included
    from the taxon names database, and only the genus and species will be searched in the sequence
    records, effectively ignoring subspecies labeling while still capturing the record. For all
    valid taxa names found, the taxon names are compiled and for each taxon all records are found.

    The sequences available for a given taxon are sorted by length (longest first). A sequence is 
    selected based on the -f flag (options: translate,length), and in the case of 
    translate the first sequence in the list that passes a translation test is selected. 
    If no sequences pass translation the first (longest) sequence is selected. 
    The translation table should be specified with the --table flag. All NCBI translation table options 
    are available, and can be selected using integers or the shortcut terms provided. If the --table 
    flag is omitted, the default will be to use the Standard translation table. 
    
    For the length strategy (-f length) the first (longest) sequence is selected. 
    For both filtering strategies, all sequences must be longer than the base pair length requirement supplied 
    (-l). If a sequence is not long enough, the search continues down the list of available 
    sequences. If no sequences are long enough, NO sequence is selected for that taxon. 

    Both strategies (translate, length) provide an objective and repeatable method to reconstruct 
    an identical supermatrix every time. If the --randomize flag is provided, the list of 
    available sequences for a given taxon is shuffled randomly, rather than sorted by length. 
    If this is used in conjunction with -f translate, a translation test and length test 
    are still performed to select a sequence. If this is used with -f length, then the sequence 
    selection is truly randomized (but the length filter must be passed). If there are 
    multiple sequences available for many taxa and loci, the -f length plus --randomize 
    option will provide a nearly infinite number of permutations of the final supermatrix 
    constructed.

    If all sequences passing the filters are desired instead of one representative sequence per
    taxon, the --allseqs flag can be used. This is useful for filtering and retaining intraspecific 
    sequence data, but it will not allow loci to be concatenated into a supermatrix downstream.

    Empirical fasta file should be labeled as 'NAME.fasta' or 'NAME.fa', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.
    		
    DEPENDENCIES: Python: BioPython.
	-----------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the fasta files to filter by species and sequence. Follow labeling format: NAME.fasta")
    parser.add_argument("-f", "--filtering", required=True, choices=["translate","length"], help="REQUIRED: Strategy for filtering sequence data.")
    parser.add_argument("-l", "--length", required=True,help="REQUIRED: Minimum number of base pairs required to keep a sequence. (an integer, ex. 150)")
    parser.add_argument("-t", "--taxa", required=True, help="REQUIRED: The full path to a text file containing all taxon names to cross-reference in the fasta file.")
    parser.add_argument("--no_subspecies", required=False, action='store_true', help="OPTIONAL: Ignore any subspecies labels in both the name database and record searches (only search binomial names).")
    parser.add_argument("--randomize", required=False, action='store_true', help="OPTIONAL: For taxa with multiple sequences, shuffle order randomly. Overrides sorting by length before searching.")
    parser.add_argument("--allseqs", required=False, action='store_true', help="OPTIONAL: For taxa with multiple sequences, select all sequences passing the filters instead of a single representative sequence.")    
    parser.add_argument("--table", default="standard", choices=["standard","vertmtdna","invertmtdna","yeastmtdna","plastid","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"], help="REQUIRED for 'translate' filter: Specifies translation table.")
    return parser.parse_args()

def split_name(string, index, delimiter):
    index = int(index)
    name = string.split(delimiter)[index]
    return name

def make_nested_dir(in_dir, string):
    dirname = in_dir+'/'+string
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    return dirname

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

def find_species(f, species, subspecies, no_subspecies):
    tax_match = set()
    with open(f, 'r') as fh_f:
        for line in fh_f:
            if line.startswith('>'):
                if '|' in line:
                    line.replace("|"," ")
                line = line.upper().strip()
                taxon_sp, taxon_ssp = get_taxon(line)
                if no_subspecies is False:
                    if taxon_match_sp(taxon_sp,species) is True:
                        if taxon_match_ssp(taxon_ssp,subspecies) is True:
                            tax_match.add(taxon_ssp.capitalize())
                        else:
                            tax_match.add(taxon_sp.capitalize())
                    elif taxon_match_sp(taxon_sp,species) is False:
                        if taxon_match_ssp(taxon_ssp,subspecies) is True:
                            tax_match.add(taxon_ssp.capitalize())
                        else:
                            pass
                #Second set of statements are when subspecies are desired
                if no_subspecies is True:
                    if taxon_match_sp(taxon_sp,species) is True:
                        tax_match.add(taxon_sp.capitalize())
                    elif taxon_match_sp(taxon_sp,species) is False:
                        pass
    species_list = sorted(tax_match)
    return species_list

def translation_val(tcode):
    tdict = {"standard":"Standard", "vertmtdna":"Vertebrate Mitochondrial", "invertmtdna":"Invertebrate Mitochondrial", "yeastmtdna":"Yeast Mitochondrial", "plastid":"11", "1":1, "2":2, "3":3, "4":4, "5":5, "6":6, "7":7, "8":8, "9":9, "10":10, "11":11, "12":12, "13":13, "14":14, "15":15, "16":16, "17":17, "18":18, "19":19, "20":20, "21":21, "22":22, "23":23, "24":24, "25":25, "26":26, "27":27, "28":28, "29":29, "30":30, "31":31}
    tval = tdict[tcode]
    return tval

def check_translation(newseq, tcode):
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
    #translate to all forward frames, constantly adjusting to keep multiples of 3
    frame1 = newseq[:].translate(table=tcode)
    f2 = newseq+n
    frame2 = f2[1:].translate(table=tcode)
    f3 = f2+n
    frame3 = f3[2:].translate(table=tcode)
    #reverse complement and translate to all frames, constantly adjusting to keep multiples of 3
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
    frames_list = [[f1stops,frame1,"Frame 1"],
                  [f2stops,frame2,"Frame 2"],
                  [f3stops,frame3,"Frame 3"],
                  [rcf1stops,rcframe1,"RC Frame 1"],
                  [rcf2stops,rcframe2,"RC Frame 2"],
                  [rcf3stops,rcframe3,"RC Frame 3"]]
        
    #create multi-conditional statement
    #first, we search all frames to see if they have 0 stop codons, if so we use that sequence frame
    #second, we search all frames to see if a single stop codon occurs in the final two positions, if so we count as pass
    match = False
    #sublist will have following contents [stop codon count, sequence in frame, frame label]
    for sublist in frames_list:
        if match is False:        
            #check if frame has no stop codons
            if sublist[0] == int(0):
                print "\t\t\t\t\t\t{0}".format(sublist[2])
                match = True

    if match is False:
        for sublist in frames_list:
            if match is False:
                #check if frame has one stop codon and it is in the last three positions of the sequence
                if sublist[0] == int(1) and '*' in sublist[1][-4:]:
                    print "\t\t\t\t\t\t{0} with single tailing stop codon".format(sublist[2])
                    match = True

    #action if all else fails
    if match is False:
        print '\t\t\t\t\t\t***Failed translation'
    return match

def species_parser_coding(f, tcode, out_dir, bp_min, randomize, species, subspecies, no_subspecies, allseqs):
    print "\n\n--------------------------------------------------------------------------------------"
    print "Processing {}".format(f)
    fasta_dict = SeqIO.index(f, "fasta")

    prefix = split_name(f,0,'.')

    if allseqs is False:
        #create the output fasta file with the a single filtered sequence per taxon
        out_fasta = "{}_single_taxon.fasta".format(prefix)
    if allseqs is True:
        #create the output fasta file with all filtered sequences for taxa
        out_fasta = "{}_all_seqs.fasta".format(prefix)        

    #create a log file containing the following headers:
    #Taxon  Accession  SeqLength   PassedTranslation  SeqsAvailable
    out_log = "{}_species_log.txt".format(prefix)
    with open(out_log, 'a') as fh_out_log:
        fh_out_log.write("Taxon\tAccession\tSeqLength\tPassedTranslation\tSeqsAvailable\n")

    #create a file with each species and all corresponding accns on a tab-delimited line 
    out_seqs = "{}_accession_list_by_species.txt".format(prefix)
    
    #create batch entrez style file which simply contains all accession numbers from the fasta file
    out_list = "{}_accession_list_for_Batch_Entrez.txt".format(prefix)

    species_list = find_species(f, species, subspecies, no_subspecies)

    #get accessions for each species in list
    for species in species_list:
        accession_list = []
        print "\t{}:".format(species)
        with open(f, 'r') as fh_f:
            for line in fh_f:
                if line.startswith('>'):
                    if '|' in line:
                        line.replace("|"," ")
                    #get accn number
                    accession = [l.strip('>') for l in line.split()][0]
                    #get taxon
                    taxon_sp, taxon_ssp = get_taxon(line)
                    #check if record matches taxon, if so, take accn number and seq length
                    #add in filter for minimum sequence length determined by user
                    if no_subspecies is False:
                        if taxon_ssp.capitalize() in species_list:
                            if taxon_ssp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                                accession_list.append([accession, len(fasta_dict[accession].seq)])
                        elif taxon_ssp.capitalize() not in species_list:
                            if taxon_sp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                                accession_list.append([accession, len(fasta_dict[accession].seq)])
                    elif no_subspecies is True:
                        if taxon_sp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                            accession_list.append([accession, len(fasta_dict[accession].seq)])
                        
        #if some sequences passed bp filter, the list will be populated
        if accession_list:
            #sort list of sublists [accn, length] by reverse length to start
            #searches with longest sequence available for the taxon
            if randomize is False:
                accession_list.sort(key=operator.itemgetter(1), reverse=True)
            #if randomize option is selected, randomly shuffle list before search
            #this screws up the sort by length feature
            elif randomize is True:
                random.shuffle(accession_list)
                
            print "\t\tSequences found: [Accession, Sequence Length]"
            for a in accession_list:
                print "\t\t\t{}".format(a)
            with open(out_seqs, 'a') as fh_out_seqs:
                fh_out_seqs.write("\n{}\t".format(species))

            if allseqs is False:
                #a way to track if a sequence has been selected yet
                kept = int(0)

                #iterate over accessions 
                for a in accession_list:
                    acc = a[0]
                    seq_len = a[1]

                    #write every accession number to batch entrez file
                    with open(out_list, 'a') as fh_out_list:
                        fh_out_list.write("{}\n".format(acc))
                    #write every accession to the species line in other log file
                    with open(out_seqs, 'a') as fh_out_seqs:
                        fh_out_seqs.write("{}\t".format(acc))

                    #if a translatable sequence hasn't been found yet
                    if kept == int(0):
                        print "\t\t\tChecking translation for {}:".format(acc)
                        #get the sequence
                        newseq = fasta_dict[acc].seq
                        #run the translation function on the sequence, will return True if passed, False if failed translation
                        translated = check_translation(newseq, tcode)

                        #if passed translation write seq to output fasta and record info in log file
                        if translated is True:
                            kept += 1
                            with open(out_fasta, 'a') as fh_out_fasta:
                                fh_out_fasta.write((fasta_dict[acc]).format("fasta"))
                            with open(out_log, 'a') as fh_out_log:
                                fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, acc, seq_len, "Y", len(accession_list)))
                        else:
                            pass
                    else:
                        pass

                #if all seqs fail translation take the first entry (the longest) and move on
                if kept == int(0):
                    print '\t\t\t\t\tNo sequences passed, taking first (longest) entry.'
                    with open(out_fasta, 'a') as fh_out_fasta:
                        fh_out_fasta.write(fasta_dict[accession_list[0][0]].format("fasta"))
                    with open(out_log, 'a') as fh_out_log:
                        fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, accession_list[0][0], accession_list[0][1], "N", len(accession_list)))

            elif allseqs is True:
                #a way to track if a sequence has been selected yet
                kept = int(0)

                #iterate over accessions 
                for a in accession_list:
                    acc = a[0]
                    seq_len = a[1]

                    #write every accession number to batch entrez file
                    with open(out_list, 'a') as fh_out_list:
                        fh_out_list.write("{}\n".format(acc))
                    #write every accession to the species line in other log file
                    with open(out_seqs, 'a') as fh_out_seqs:
                        fh_out_seqs.write("{}\t".format(acc))

                    print "\t\t\tChecking translation for {}:".format(acc)
                    #get the sequence
                    newseq = fasta_dict[acc].seq
                    #run the translation function on the sequence, will return True if passed, False if failed translation
                    translated = check_translation(newseq, tcode)

                    #if passed translation write seq to output fasta and record info in log file
                    if translated is True:
                        kept += 1
                        with open(out_fasta, 'a') as fh_out_fasta:
                            fh_out_fasta.write((fasta_dict[acc]).format("fasta"))
                        with open(out_log, 'a') as fh_out_log:
                            fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, acc, seq_len, "Y", len(accession_list)))
                    else:
                        pass

                #if all seqs fail translation take the first entry (the longest) and move on
                if kept == int(0):
                    print '\t\t\t\t\tNo sequences passed, taking first (longest) entry.'
                    with open(out_fasta, 'a') as fh_out_fasta:
                        fh_out_fasta.write(fasta_dict[accession_list[0][0]].format("fasta"))
                    with open(out_log, 'a') as fh_out_log:
                        fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, accession_list[0][0], accession_list[0][1], "N", len(accession_list)))
                
        #if list is empty no seqs passed bp filter
        elif not accession_list:
            print "\t\t***No sequences passed the {} bp minimum filter***\n".format(bp_min)
                
    shutil.move(out_fasta, out_dir)
    shutil.move(out_log, out_dir)
    shutil.move(out_seqs, out_dir)
    shutil.move(out_list, out_dir)

def species_parser_noncoding(f, out_dir, bp_min, randomize, species, subspecies, no_subspecies, allseqs):
    print "\n\n--------------------------------------------------------------------------------------"
    print "Processing {}".format(f)
    fasta_dict = SeqIO.index(f, "fasta")

    prefix = split_name(f,0,'.')

    if allseqs is False:
        #create the output fasta file with the a single filtered sequence per taxon
        out_fasta = "{}_single_taxon.fasta".format(prefix)
    if allseqs is True:
        #create the output fasta file with all filtered sequences for taxa
        out_fasta = "{}_all_seqs.fasta".format(prefix)        

    #create a log file containing the following headers:
    #Taxon  Accession  SeqLength   PassedTranslation  SeqsAvailable
    out_log = "{}_species_log.txt".format(prefix)
    with open(out_log, 'a') as fh_out_log:
        fh_out_log.write("Taxon\tAccession\tSeqLength\tPassedTranslation\tSeqsAvailable\n")

    #create a file with two headers, species and accn
    out_seqs = "{}_accession_list_by_species.txt".format(prefix)

    #create batch entrez style file which simply contains all accession numbers from the fasta file
    out_list = "{}_accession_list_for_Batch_Entrez.txt".format(prefix)

    species_list = find_species(f, species, subspecies, no_subspecies)
    
    #get accessions for each species in list
    for species in species_list:
        accession_list = []
        print "\t{}:".format(species)
        with open(f, 'r') as fh_f:
            for line in fh_f:
                if line.startswith('>'):
                    if '|' in line:
                        line.replace("|"," ")
                    #get accn number
                    accession = [l.strip('>') for l in line.split()][0]
                    #get taxon
                    taxon_sp, taxon_ssp = get_taxon(line)
                    #check if record matches taxon, if so, take accn number and seq length
                    #add in filter for minimum sequence length determined by user
                    if no_subspecies is False:
                        if taxon_ssp.capitalize() in species_list:
                            if taxon_ssp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                                accession_list.append([accession, len(fasta_dict[accession].seq)])
                        elif taxon_ssp.capitalize() not in species_list:
                            if taxon_sp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                                accession_list.append([accession, len(fasta_dict[accession].seq)])
                    elif no_subspecies is True:
                        if taxon_sp.capitalize() == species and len(fasta_dict[accession].seq) >= int(bp_min):
                            accession_list.append([accession, len(fasta_dict[accession].seq)])

        #if some sequences passed bp filter, the list will be populated
        if accession_list:
            with open(out_seqs, 'a') as fh_out_seqs:
                fh_out_seqs.write("\n{}\t".format(species))
            
            #sort list of sublists [accn, length] by reverse length        
            if randomize is False:
                accession_list.sort(key=operator.itemgetter(1), reverse=True)
            #if randomize option is selected, randomly shuffle list before search
            #this screws up the sort by length feature
            elif randomize is True:
                random.shuffle(accession_list)

            print "\t\tSequences found: [Accession, Sequence Length]"
            for a in accession_list:
                print "\t\t\t{}".format(a)
                #write each accession number to batch entrez file
                with open(out_list, 'a') as fh_out_list:
                    fh_out_list.write("{}\n".format(a[0]))
                #write every accession to the species line in other log file
                with open(out_seqs, 'a') as fh_out_seqs:
                    fh_out_seqs.write("{}\t".format(a[0]))
                if allseqs is True:
                    with open(out_fasta, 'a') as fh_out_fasta:
                        fh_out_fasta.write(fasta_dict[a[0]].format("fasta"))
                    with open(out_log, 'a') as fh_out_log:
                        fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, a[0], a[1], "NA", len(accession_list)))
                    
            if allseqs is False:
                #this is easy, just take first (longest) sequence in list of sublists [accn, seqlength]
                with open(out_fasta, 'a') as fh_out_fasta:
                    fh_out_fasta.write(fasta_dict[accession_list[0][0]].format("fasta"))
                with open(out_log, 'a') as fh_out_log:
                    fh_out_log.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(species, accession_list[0][0], accession_list[0][1], "NA", len(accession_list)))
                                
        #if list is empty no seqs passed bp filter
        elif not accession_list:
            print "\t\t***No sequences passed the {} bp minimum filter***\n".format(bp_min)
                        
    shutil.move(out_fasta, out_dir)
    shutil.move(out_log, out_dir)
    shutil.move(out_seqs, out_dir)
    shutil.move(out_list, out_dir)

#-----------------------------------------------------------------------------------------

def main():
    args = get_args()
    species, subspecies = parse_taxa(args.taxa, args.no_subspecies)
    os.chdir(args.in_dir)
    fasta_list = sorted([f for f in os.listdir('.') if f.endswith('.fasta') or f.endswith('.fa')])
    print "Found {} fasta files to filter...\n\n".format(len(fasta_list))
    filter_dir = make_nested_dir(args.in_dir, 'Species_Filtering_Results')

    if args.filtering == "translate":
        tb = datetime.now()
        if args.table is None:
            tcode = "Standard"
        else:
            tcode = translation_val(args.table)
        for f in fasta_list:
            tbf = datetime.now()
            species_parser_coding(f, tcode, filter_dir, args.length, args.randomize, species, subspecies, args.no_subspecies, args.allseqs)
            tff = datetime.now()
            tef = tff - tbf
            print "\n\nTotal time to process {0}: {1} (H:M:S)".format(f,tef)
            print "--------------------------------------------------------------------------------------\n\n"
        tf = datetime.now()
        te = tf - tb
        print "\n\nTotal time to filter all files: {1} (H:M:S)\n\n".format(f,te)

    elif args.filtering == "length":
        tb = datetime.now()
        for f in fasta_list:
            tbf = datetime.now()
            species_parser_noncoding(f, filter_dir, args.length, args.randomize, species, subspecies, args.no_subspecies, args.allseqs)
            tff = datetime.now()
            tef = tff - tbf
            print "\n\nTotal time to process {0}: {1} (H:M:S)".format(f,tef)
            print "--------------------------------------------------------------------------------------\n\n"
        tf = datetime.now()
        te = tf - tb
        print "\n\nTotal time to filter all files: {0} (H:M:S)\n\n".format(te)
            
if __name__ == '__main__':
    main()
