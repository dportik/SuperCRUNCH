''''
SuperCRUNCH: Relabel_Fasta module

Usage: Usage: python Relabel_Fasta.py -i [directory with all fasta files] (REQUIRED)
                                      -r [strategy for relabeling: "species", "accession", "species_acc"] REQUIRED
                                      -s [full path to text file with subspecies names to include] (OPTIONAL)

    Relabel_Fasta: Process fasta files to relabel record using one of three strategies:
    species, accession, species_acc. With 'species', the record description is split by 
    spaces, exclude the accession number, then join the second and third elements of 
    the split line with an underscore. This generally corresponds to the genus and 
    species if records are labeled properly. With 'accession', the accession number
    will be used to label the sequence record. With 'species_acc', the record is 
    labeled by taxon and accession number. In all strategies, any spaces are replaced
    by underscores. If the optional -s (--subspecies) flag is included with a text
    file containing subspecies (trinomial) names, then the species component of each
    of the above options will include the trinomial if appropriate.

    Example input fasta file:

    >JN881132.1 Daboia russelii activity-dependent neuroprotector (ADNP) gene, partial cds
    ...Sequence....
    >KU765220.1 Sceloporus undulatus voucher ADL182 activity-dependent neuroprotector (adnp) gene, partial cds
    ...Sequence....

    Output using -r 'species':

    >Daboia_russelii
    ...Sequence....
    >Sceloporus_undulatus
    ...Sequence....

    Output using -r 'accession':

    >JN881132.1
    ...Sequence....
    >KU765220.1
    ...Sequence....

    Output using -r 'species_acc':

    >Daboia_russelii_JN881132.1
    ...Sequence....
    >Sceloporus_undulatus_KU765220.1
    ...Sequence....

    Empirical fasta file should be labeled as 'NAME.fasta', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

[-s] The input taxon name file should simply contain a list of taxon names, one on each line. 
    Although this file is only necessary for getting subspecies taxonomic names, it can also 
    contain binomial names (genus and species) that will simply be ignored. For all names the 
    genus, species, and subspecies should be separated by a space. All taxon names are converted 
    to uppercase for searching, so names are NOT case-sensitive.
 
	Example of file structure for taxon information (showing species and subspecies examples):
	
    Varanus acanthurus
    Varanus albigularis albigularis
    Varanus albigularis microstictus
    Varanus auffenbergi
    Varanus bangonorum
    Varanus baritji
    Varanus beccarii
    
    Here is how the subspecies flag would affect the following labeling of an example input fasta file:

    >JN881132.1 Varanus acanthurus activity-dependent neuroprotector (ADNP) gene, partial cds
    ...Sequence....
    >KU765220.1 Varanus albigularis albigularis voucher ADL182 activity-dependent neuroprotector (adnp) gene, partial cds
    ...Sequence....

    Output using -r 'species':

    >Varanus_acanthurus
    ...Sequence....
    >Varanus_albigularis_albigularis
    ...Sequence....

    Output using -r 'accession':

    >JN881132.1
    ...Sequence....
    >KU765220.1
    ...Sequence....

    Output using -r 'species_acc':

    >Varanus_acanthurus_JN881132.1
    ...Sequence....
    >Varanus_albigularis_albigularis_KU765220.1
    ...Sequence....


    Output files:
    
        [fasta name]_relabeled.fasta - Contains all the relabeled sequences from original input fasta.
        								No modifications have been made to the actual sequences.

-------------------------
For Python 2.7
Dependencies: None
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
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Relabel_Fasta: Process fasta files to relabel record using one of three strategies:
    species, accession, species_acc. With 'species', the record description is split by 
    spaces, exclude the accession number, then join the second and third elements of 
    the split line with an underscore. This generally corresponds to the genus and 
    species if records are labeled properly. With 'accession', the accession number
    will be used to label the sequence record. With 'species_acc', the record is 
    labeled by taxon and accession number. In all strategies, any spaces are replaced
    by underscores. If the optional -s (--subspecies) flag is included with a text
    file containing subspecies (trinomial) names, then the species component of each
    of the above options will include the trinomial if appropriate.

    Empirical fasta file should be labeled as 'NAME.fasta', where NAME represents the
    gene/locus. The NAME portion should not contain any periods or spaces, but can contain
    underscores. Output files are labeled using a prefix identical to NAME.

    Output files:
	[NAME]_relabeled.fasta - Contains relabeled description lines for all records from original input fasta.
	DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input fasta files. Follow labeling format: NAME.fasta")
    parser.add_argument("-r", "--relabel", required=True, choices=["species", "accession", "species_acc"], help="REQUIRED: The strategy for relabeling seequence records.")
    parser.add_argument("-s", "--subspecies", required=False, default=None, help="OPTIONAL: The full path to a text file containing all subspecies names to cross-reference in the fasta file.")
    return parser.parse_args()


def parse_taxa(f):
    print "\nParsing taxon information from {}.".format(f)
    with open(f, 'r') as fh_f:
        subspecies_set = set([line.upper().strip().replace(" ","_") for line in fh_f if len(line.split()) == int(3)])
    subspecies = sorted(subspecies_set)
    print "\tFound {} subspecies names to include for relabeling.\n".format(len(subspecies))
    return subspecies

def get_taxon(line):
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split() if line.split() >= int(3)][1:3]
    taxon_sp = "_".join(parts1)
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split() if line.split() >= int(4)][1:4]
    taxon_ssp = "_".join(parts2)
    return taxon_sp, taxon_ssp

def relabel_species(f,out_dir,subspecies):
    print "\nRelabeling sequence records with taxon name in {}.".format(f)
    out_f = "{}_relabeled.fasta".format(f.split('.')[0])
    with open(out_f, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):
                    if '|' in line:
                        line.replace("|"," ")
                    taxon_sp, taxon_ssp = get_taxon(line)
                    if subspecies is not None:
                        if taxon_ssp in subspecies:
                            newline = ">{}\n".format(taxon_ssp.capitalize())
                        else:
                            newline = ">{}\n".format(taxon_sp.capitalize())
                    else:
                        newline = ">{}\n".format(taxon_sp.capitalize())
                    fh_out.write(newline)
                else:
                    fh_out.write(line+'\n')
    shutil.move(out_f, out_dir)
    print "\tDone."
    
def relabel_accession(f,out_dir):
    print "\nRelabeling sequence records with accession number in {}.".format(f)
    out_f = "{}_relabeled.fasta".format(f.split('.')[0])
    with open(out_f, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):
                    if '|' in line:
                        line.replace("|"," ")
                    acc = [l.strip('>') for l in line.split()][0]
                    newline = ">{}\n".format(acc)
                    fh_out.write(newline)
                else:
                    fh_out.write(line+'\n')
    shutil.move(out_f, out_dir)
    print "\tDone."

def relabel_species_acc(f,out_dir,subspecies):
    print "\nRelabeling sequence records with taxon name and accession number in {}.".format(f)
    out_f = "{}_relabeled.fasta".format(f.split('.')[0])
    with open(out_f, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):
                    if '|' in line:
                        line.replace("|"," ")
                    acc = [l.strip('>') for l in line.split()][0]
                    taxon_sp, taxon_ssp = get_taxon(line)
                    if subspecies is not None:
                        if taxon_ssp in subspecies:
                            newline = ">{0}_{1}\n".format(taxon_ssp.capitalize(),acc)
                        else:
                            newline = ">{0}_{1}\n".format(taxon_sp.capitalize(),acc)
                    else:
                        newline = ">{0}_{1}\n".format(taxon_sp.capitalize(),acc)
                    fh_out.write(newline)
                else:
                    fh_out.write(line+'\n')
    shutil.move(out_f, out_dir)
    print "\tDone."
    
def main():
    tb = datetime.now()
    args = get_args()
    os.chdir(args.in_dir)
    if args.relabel == "species":
        out_dir = 'Relabeled_Fasta_Files_Species'
    elif args.relabel == "accession":
        out_dir = 'Relabeled_Fasta_Files_Accession'
    elif args.relabel == "species_acc":
        out_dir = 'Relabeled_Fasta_Files_SpeciesAccession'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    f_list = sorted([f for f in os.listdir('.') if f.endswith(".fasta") or f.endswith(".fa")])

    subspecies = None
    if args.subspecies is not None:
        subspecies = parse_taxa(args.subspecies)
    
    for f in f_list:
        if args.relabel == "species":
            relabel_species(f,out_dir,subspecies)
        elif args.relabel == "accession":
            relabel_accession(f,out_dir)
        elif args.relabel == "species_acc":
            relabel_species_acc(f,out_dir,subspecies)
            
    tf = datetime.now()
    te = tf - tb
    print "\n\n\nFinished relabeling alignments. Elapsed time: {} (H:M:S)\n\n".format(te)

if __name__ == '__main__':
    main()
