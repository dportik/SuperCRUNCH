'''
SuperCRUNCH: Make_Acc_Table module

    Make_Acc_Table: A tool for extracting taxon and accession number 
    information from a directory of fasta files to create a table
    of GenBank accession numbers.

    A comprehensive list of taxa is collected across all fasta files 
    in the directory. The taxon names recovered are two-part by default 
    (e.g., Genus species). To include three-part names (e.g., Genus 
    species subspecies), the optional -s flag can be used. This requires
    the full path to a taxon names file. This file can contain a mix of 
    two-part and three-part names, and can be the same file used for 
    earlier steps, but only the three-part names are used to aid the 
    searches. If this is a vouchered data set (produced from the 
    Filter_Seqs_and_Species module using the voucherize option), the
    --voucherize flag should be used here to construct the taxon labels 
    with the voucher information present. This will produce an accession 
    table where each entry is a taxon + specific sample, rather than just 
    a taxon label. If a vouchered data set is run without the --voucherize 
    flag, an error will be thrown indicating duplicate taxon labels are 
    present within a single fasta file.
    
    The final table is written with taxon labels in the first column 
    (sorted alphabetically) and the remaining columns representing the 
    fasta files (also sorted alphabetically), with rows filled with 
    accession numbers or dashes depending on the presence or absence 
    of that taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files, as long as
    they contain full description lines for sequences.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

-------------------------
Compatible with Python 2.7 & 3.7
Dependencies: 
    None
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
January 2019
Distributed under the 
GNU General Public Lincense
'''
import os
import shutil
import argparse
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Make_Acc_Table: A tool for extracting taxon and accession number 
    information from a directory of fasta files to create a table
    of GenBank accession numbers. Can be used for supermatrix datasets 
    and vouchered population-level datasets. 

    A comprehensive list of taxa is collected across all fasta files 
    in the directory. The taxon names recovered are two-part by default 
    (e.g., Genus species). To include three-part names (e.g., Genus 
    species subspecies), the optional -s flag can be used. This requires
    the full path to a taxon names file. This file can contain a mix of 
    two-part and three-part names, and can be the same file used for 
    earlier steps, but only the three-part names are used to aid the 
    searches. If this is a vouchered data set (produced from the 
    Filter_Seqs_and_Species module using the voucherize option), the
    --voucherize flag should be used here to construct the taxon labels 
    with the voucher information that is present. This will produce an accession 
    table where each entry is a taxon + specific sample, rather than just 
    a taxon label. If a vouchered data set is run without the --voucherize 
    flag, an error will be thrown indicating duplicate taxon labels are 
    present within a single fasta file.
    
    The final table is written with taxon labels in the first column 
    (sorted alphabetically) and the remaining columns representing the 
    fasta files (also sorted alphabetically), with rows filled with 
    accession numbers or dashes depending on the presence or absence 
    of that taxon for a given locus.

    This tool can be used on unaligned or aligned fasta files, as long as
    they contain full description lines for sequences.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.
    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the input fasta files.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-s", "--subspecies",
                            required=False,
                            default=None,
                            help="OPTIONAL: If subspecies names are to be included, this "
                            "should be the full path to a text file containing all "
                            "subspecies names to cross-reference in the fasta file. "
                            "This can be the same taxon list used in earlier steps.")

    parser.add_argument("--voucherize",
                            required=False,
                            action='store_true',
                            help="OPTIONAL: If a 'Voucher_[ID]' entry is included in the "
                            "sequence description, append the ID component to the end of "
                            "the species/subspecies label (e.g., Hyperolius_nitidulus_MVZ236473). "
                            "The Voucher_[ID] field is generated by the Parse_Loci.py module.")
    
    return parser.parse_args()        

def parse_taxa(f):
    """
    Retrieve subspecies (three-part) names ONLY from user supplied taxon 
    names file (f). Names are joined by underscores and converted to
    capital case. For example:
    Original name - Hyperolius balfouri viridistriatus
    Processed name - Hyperolius_balfouri_viridistriatus
    Returns a sorted list of these names.
    """
    print("\nGetting subspecies from {}.".format(f.split('/')[-1]))
    
    with open(f, 'r') as fh:
        subspecies_set = set([line.upper().strip().replace(" ","_").capitalize() for line in fh
                                  if len(line.split()) == int(3)])
        
    subspecies = sorted(subspecies_set)
    print("\tFound {} subspecies names to include.".format(len(subspecies)))
    
    return subspecies

def get_taxon(line):
    """
    Retrieve the taxon name from '>' line in a fasta file.
    Will fetch the species (two-part) and subspecies (three-part)
    names. If voucherize argument was supplied, will attempt to 
    identify a 'Voucher_SOMETHING' element in the description line
    and add the 'SOMETHING' component to the end of the species
    and subspecies names. Names are joined by underscores and 
    converted to capital case. 
    """
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if len(line.split()) >= int(3)][1:3]
    taxon_sp = ("_".join(parts1)).capitalize()
    
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if len(line.split()) >= int(4)][1:4]
    taxon_ssp = ("_".join(parts2)).capitalize()

    taxon_spv = ""
    taxon_sspv = ""
    
    for part in line.split():
        if part.startswith('Voucher_'):
            taxon_spv = "{}_{}".format(taxon_sp, part.replace('Voucher_',''))
            taxon_sspv = "{}_{}".format(taxon_ssp, part.replace('Voucher_',''))
            
    if not taxon_spv:
        taxon_spv = taxon_sp
        
    if not taxon_sspv:
        taxon_sspv = taxon_ssp
    
    return taxon_sp, taxon_ssp, taxon_spv, taxon_sspv

def get_taxa_fasta(f, subspecies_list, voucherize):
    """
    Retrieve the names following the > in all lines
    from a fasta file (f) based on subspecies option
    and voucherize option, add names to the set,
    then return the set.
    """
    taxa_set = set()
    
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith(">"):
                taxon_sp, taxon_ssp, taxon_spv, taxon_sspv = get_taxon(line)

                if voucherize:
                    if subspecies_list:
                        if taxon_ssp in subspecies_list:
                            taxa_set.add(taxon_sspv)
                        else:
                            taxa_set.add(taxon_spv)
                    else:
                        taxa_set.add(taxon_spv)

                else:
                    if subspecies_list:
                        if taxon_ssp in subspecies_list:
                            taxa_set.add(taxon_ssp)
                        else:
                            taxa_set.add(taxon_sp)
                    else:
                        taxa_set.add(taxon_sp)
    
    return taxa_set
    
def collect_taxa(path_prefix_list, subspecies_list, voucherize):
    """
    Collect taxon names from all alignment files (flist)
    based on subspecies option. Create a set of names 
    from each file, add to a larger set, then return 
    the larger set as a sorted list.
    """
    taxon_set = set()
    
    for item in path_prefix_list:
        fset = get_taxa_fasta(item[0], subspecies_list, voucherize)
        taxon_set.update(fset)
        
    taxon_list = sorted(taxon_set)
    
    return taxon_list


def dict_taxon_accs(f, subspecies_list, voucherize):
    """
    Function to convert fasta file (f) into
    dictionary structure with taxon label as key
    and the accession number as the value, based 
    on the subspecies and voucherize options.
    """
    fdict = {}
    fname = f.split('/')[-1]
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line.split()[0].strip('>')
                taxon_sp, taxon_ssp, taxon_spv, taxon_sspv = get_taxon(line)

                if voucherize:
                    if subspecies_list:
                        if taxon_ssp in subspecies_list:
                            fdict[taxon_sspv] = acc
                        else:
                            fdict[taxon_spv] = acc
                    else:
                        fdict[taxon_spv] = acc

                else:
                    if subspecies_list:
                        if taxon_ssp in subspecies_list:
                            if taxon_ssp in fdict:
                                raise ValueError("\n\n\n***Warning: Multiple accession numbers detected for"
                                                     " '{}' in file {}.\n\tTry using the --voucherize option if this is a "
                                                     "vouchered data set, or check your files for errors.\n\n".format(taxon_ssp, fname))
                            else:
                                fdict[taxon_ssp] = acc
                        else:
                            if taxon_sp in fdict:
                                raise ValueError("\n\n\n***Warning: Multiple accession numbers detected for"
                                                     " '{}' in file {}.\n\tTry using the --voucherize option if this is a "
                                                     "vouchered data set, or check your files for errors.\n\n".format(taxon_sp, fname))
                            else:
                                fdict[taxon_sp] = acc
                    else:
                        if taxon_sp in fdict:
                            raise ValueError("\n\n\n***Warning: Multiple accession numbers detected for"
                                                " '{}' in file {}.\n\tTry using the --voucherize option if this is a "
                                                "vouchered data set, add subspecies using -s if this contains subspecies, "
                                                 "or check your files for errors.\n\n".format(taxon_sp, fname))
                        else:
                            fdict[taxon_sp] = acc

    return fdict

def taxa_acc_dict(dict_list, taxa_list):
    """
    Create a dictionary where keys are taxa from
    taxa_list and values are accession numbers or
    missing data values (-) from all the fasta files,
    joined by tabs. Because the dict_list
    is sorted in the same order as the fasta files,
    the accession numbers are written in the same
    ordered sequence.
    """
    #initiate empty dictionary
    writing_dict = {}

    #iterate over all taxa in taxon list
    for taxon in taxa_list:
        #create key with empty string placeholder value
        writing_dict[taxon] = ""
        #iterate over dictionary list
        for d in dict_list:
            #check if current taxon name is a key in this
            #particular dictionary (d)
            if taxon in d:
                #if so, write the value (accession) of this
                #dictionary (d) as part of the new value in 
                #writing_dict for this taxon
                writing_dict[taxon] += "{}\t".format(d[taxon])
            #if current taxon name is not a key in this dictionary (d)
            else:
                #write missing data symbol as part of the new
                #value in writing_dict for this taxon                
                writing_dict[taxon] += "{}\t".format("-")
                
    return writing_dict

def write_acc_table(taxa_list, writing_dict, path_prefix_list):
    """
    Writes an accession table for all entries of the taxon list.
    Uses the dictionary structure (writing_dict) in which keys are taxa
    and values are the previously concatenated accession numbers. 
    """
    with open("GenBank_Accession_Table.txt", 'a') as fh:
        #create column labels first
        #write a column header called Taxon
        fh.write('Taxon\t')
        #then create a column header for each fasta file name
        for item in path_prefix_list:
            fh.write('{}\t'.format(item[1]))
        fh.write('\n')
        #now write the entries for each taxon
        for taxon in taxa_list:
            fh.write("{0}\t{1}\n".format(taxon.replace("_"," "),
                                             writing_dict[taxon].strip('\t')))

def write_taxon_summary(taxa_list, writing_dict):
    """
    Write a simple file with two columns, the taxon label
    and the number of loci found for that taxon. 
    """
    with open("Taxon_Summary.txt", 'a') as fh:
        #write headers
        fh.write("{}\t{}\n".format("Taxon", "Loci"))
        #iterate over taxon list
        for taxon in taxa_list:
            #break concatenated accession string into list
            acc_list = writing_dict[taxon].strip('\t').split('\t')
            #create new list without missing data
            loci_present = [a for a in acc_list if a != '-']
            #write to file
            fh.write("{}\t{}\n".format(taxon, len(loci_present)))
    
def main():
    tb = datetime.now()
    args = get_args()
    
    if args.subspecies:
        print("\nUsing subspecies labels found in: {}".format(args.subspecies.split('/')[-1]))
        subspecies_list = parse_taxa(args.subspecies)
        
    else:
        subspecies_list = None
        
    if args.voucherize:
        print("\nVoucherize option selected.")
       
    os.chdir(args.indir)

    #get paths to fasta files as well as their prefix names as sublists in larger list
    path_prefix_list = sorted([[os.path.abspath(f), f.split(".")[0]] for f in os.listdir('.')
                        if f.endswith((".fasta", ".fa"))])
    print("\nFound {:,} files to use for accession table.".format(len(path_prefix_list)))

    #find all taxon names in the fasta files
    taxa_list = collect_taxa(path_prefix_list, subspecies_list, args.voucherize)
    print("\nFound {:,} unique taxon labels across fasta files.".format(len(taxa_list)))

    #find all accession numbers for each taxon label
    print("\nGathering accession numbers for all taxa.")
    dict_list = [dict_taxon_accs(f[0], subspecies_list, args.voucherize) for f in path_prefix_list]
    writing_dict = taxa_acc_dict(dict_list, taxa_list)

    #write output files
    os.chdir(args.outdir)
    print("\nWriting accession table.")
    write_acc_table(taxa_list, writing_dict, path_prefix_list)
    print("\nWriting taxon summary file.")
    write_taxon_summary(taxa_list, writing_dict)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()

