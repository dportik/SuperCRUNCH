''''
SuperCRUNCH: Relabel_Fasta module

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
    >KU765220.1 Sceloporus undulatus voucher ADL182 activity-dependent neuroprotector ...
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

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

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
    
    Here is how the subspecies flag would affect the following labeling of an example 
    input fasta file:

    >JN881132.1 Varanus acanthurus activity-dependent neuroprotector (ADNP) gene, ...
    ...Sequence....
    >KU765220.1 Varanus albigularis albigularis voucher ADL182 activity-dependent ...
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
    
        [fasta name]_relabeled.fasta - Contains all the relabeled sequences from 
                                    original input fasta. No modifications have been made 
                                    to the actual sequences.

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
import operator
import argparse
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Relabel_Fasta: Process fasta files to relabel record using one of three 
    strategies: species, accession, species_acc. With 'species', the record 
    description is split by spaces, exclude the accession number, then join 
    the second and third elements of the split line with an underscore. This 
    generally corresponds to the genus and species if records are labeled 
    properly. With 'accession', the accession number will be used to label 
    the sequence record. With 'species_acc', the record is labeled by taxon 
    and accession number. In all strategies, any spaces are replaced by 
    underscores. If the optional -s (--subspecies) flag is included with a 
    text file containing subspecies (trinomial) names, then the species 
    component of each of the above options will include the trinomial if 
    appropriate. 
    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.
    DEPENDENCIES: None.
    ---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which "
                            "contains the input fasta files. Follow labeling format:"
                            "NAME.fasta")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-r", "--relabel",
                            required=True,
                            choices=["species", "accession", "species_acc"],
                            help="REQUIRED: The strategy for relabeling seequence records.")
    
    parser.add_argument("-s", "--subspecies",
                            required=False,
                            default=None,
                            help="OPTIONAL: Allows subspecies names to be included in the "
                            "relabeling. This flag requires a full path to a text file containing "
                            "all subspecies names to cross-reference in the fasta file.")
    
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
    Retrieve taxon names from user supplied taxon names file (f).
    Here, will ONLY target subspecies names. Returns a list of 
    subspecies, with all names in uppercase.
    """
    print("\nParsing taxon information from {}.".format(f))
    
    with open(f, 'r') as fh:
        subspecies = sorted(set([line.upper().strip().replace(" ","_") for line in fh
                                  if len(line.split()) == int(3)]))
        
    print("\tFound {} subspecies names to include for relabeling.\n".format(len(subspecies)))
    
    return subspecies

def get_taxon(line):
    """
    Retrieve the taxon name from '>' line in a fasta file.
    Will fetch the species (binomial) and subspecies (trinomial)
    names. If voucherize argument was supplied, will attempt to 
    identify a 'Voucher_SOMETHING' element in the description line
    and add the 'SOMETHING' component to the end of the species
    and subspecies names.
    """
    parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if len(line.split()) >= int(3)][1:3]
    taxon_sp = ("_".join(parts1)).capitalize()
    
    parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.upper().split()
                  if len(line.split()) >= int(4)][1:4]
    
    if parts2[2].isalpha():
        taxon_ssp = ("_".join(parts2)).capitalize()
    else:
        if "_" in parts2[2]:
            substrings = parts2[2].split('_')
            if substrings[0].isalpha() == True and substrings[1].isalpha() == False:
                taxon_ssp = "_".join([parts2[0].capitalize(), parts2[1].lower(), substrings[0].lower(), substrings[1].upper()])
            else:
                taxon_ssp = "_".join([parts2[0].capitalize(), parts2[1].lower(), parts2[2].upper()])
        else:
            taxon_ssp = "_".join([parts2[0].capitalize(), parts2[1].lower(), parts2[2].upper()])

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

def relabel_species(f, fdir, subspecies, voucherize):
    """
    Relabel sequence records using a taxon label.
    """
    label_key = []
    print("\nRelabeling sequence records with taxon name in {}.".format(f))
    outname = "{}_relabeled.fasta".format(f.split('.')[0])
    
    with open(outname, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):                        
                    acc = [l.strip('>') for l in line.split()][0]
                    descrip = " ".join([l.strip() for l in line.split()][1:])
                    
                    taxon_sp, taxon_ssp, taxon_spv, taxon_sspv = get_taxon(line)

                    if voucherize is True:
                        if subspecies:
                            if taxon_ssp.upper() in subspecies:
                                newline = ">{}\n".format(taxon_sspv)
                                label_key.append([acc, taxon_ssp.replace('_',' '), descrip])
                            else:
                                newline = ">{}\n".format(taxon_spv)
                                label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        else:
                            newline = ">{}\n".format(taxon_spv)
                            label_key.append([acc, taxon_sp.replace('_',' '), descrip])
                        fh_out.write(newline)
                        
                    elif voucherize is False:
                        if subspecies:
                            if taxon_ssp.upper() in subspecies:
                                newline = ">{}\n".format(taxon_ssp)
                                label_key.append([acc, taxon_ssp.replace('_',' '), descrip])
                            else:
                                newline = ">{}\n".format(taxon_sp)
                                label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        else:
                            newline = ">{}\n".format(taxon_sp)
                            label_key.append([acc, taxon_sp.replace('_',' '), descrip])
                        fh_out.write(newline)
                    
                else:
                    fh_out.write("{}\n".format(line))
                    
    shutil.move(outname, fdir)
    print("\tDone.")
    
    return label_key
    
def relabel_accession(f, fdir, subspecies, voucherize):
    """
    Relabel sequence records using an accession number.
    """
    label_key = []
    print("\nRelabeling sequence records with accession number in {}.".format(f))
    outname = "{}_relabeled.fasta".format(f.split('.')[0])
    
    with open(outname, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):
                    acc = [l.strip('>') for l in line.split()][0]
                    descrip = " ".join([l.strip() for l in line.split()][1:])
                    
                    newline = ">{}\n".format(acc)
                    fh_out.write(newline)
                    
                    taxon_sp, taxon_ssp, taxon_spv, taxon_sspv = get_taxon(line)
                    
                    if subspecies:
                        if taxon_ssp in subspecies:
                            label_key.append([acc, taxon_ssp.replace('_',' '), descrip])
                        else:
                            label_key.append([acc, taxon_sp.replace('_',' '), descrip])
                            
                    else:
                        label_key.append([acc, taxon_sp.replace('_',' '), descrip])
                        
                else:
                    fh_out.write("{}\n".format(line))
                    
    shutil.move(outname, fdir)
    print("\tDone.")
    
    return label_key

def relabel_species_acc(f, fdir, subspecies, voucherize):
    """
    Relabel sequence records using a taxon label and accession number.
    """
    label_key = []
    print("\nRelabeling sequence records with taxon name and accession number in {}.".format(f))
    outname = "{}_relabeled.fasta".format(f.split('.')[0])
    
    with open(outname, 'a') as fh_out:
        with open(f, 'r') as fh_temp:
            for line in fh_temp:
                line = line.strip()
                if line.startswith(">"):
                    acc = [l.strip('>') for l in line.split()][0]
                    descrip = " ".join([l.strip() for l in line.split()][1:])
                    
                    taxon_sp, taxon_ssp, taxon_spv, taxon_sspv = get_taxon(line)
                    
                    if voucherize is True:
                        if subspecies:
                            if taxon_ssp in subspecies:
                                newline = ">{0}_{1}\n".format(taxon_sspv, acc)
                                label_key.append([acc, taxon_ssp.replace('_',' '), descrip])
                            else:
                                newline = ">{0}_{1}\n".format(taxon_spv, acc)
                                label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        else:
                            newline = ">{0}_{1}\n".format(taxon_spv, acc)
                            label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        fh_out.write(newline)
                        
                    elif voucherize is False:
                        if subspecies:
                            if taxon_ssp in subspecies:
                                newline = ">{0}_{1}\n".format(taxon_ssp, acc)
                                label_key.append([acc, taxon_ssp.replace('_',' '), descrip])
                            else:
                                newline = ">{0}_{1}\n".format(taxon_sp, acc)
                                label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        else:
                            newline = ">{0}_{1}\n".format(taxon_sp, acc)
                            label_key.append([acc, taxon_sp.replace('_',' '), descrip])

                        fh_out.write(newline)
                    
                else:
                    fh_out.write("{}\n".format(line))
                    
    shutil.move(outname, fdir)
    print("\tDone.")
    
    return label_key

def write_label_key(label_key, f, ldir):
    outname = "{}_label_key.txt".format(f.split('.')[0])
    
    with open(outname, 'a') as fh_out:
        fh_out.write("{}\t{}\t{}\n".format("Accession", "Taxon", "Description"))
        
    label_key.sort(key=operator.itemgetter(1))
    
    with open(outname, 'a') as fh_out:
        for l in label_key:
            fh_out.write("{}\t{}\t{}\n".format(l[0], l[1], l[2]))
            
    shutil.move(outname, ldir)


def make_dirs(outdir, relabel):
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    os.chdir(outdir)
    curpath = os.getcwd()
    
    maindir = os.path.join(curpath, "Relabeled-by-{}".format(relabel))
    if not os.path.exists(maindir):
        os.mkdir(maindir)

    ldir = os.path.join(curpath, "Relabeled-by-{}".format(relabel), "Relabeling-key-files")
    if not os.path.exists(ldir):
        os.mkdir(ldir)
        
    fdir = os.path.join(curpath, "Relabeled-by-{}".format(relabel), "Relabeled-fasta-files")
    if not os.path.exists(fdir):
        os.mkdir(fdir)
        
    return fdir, ldir

def main():
    tb = datetime.now()
    args = get_args()

    fdir, ldir = make_dirs(args.outdir, args.relabel)

    os.chdir(args.indir)
    flist = sorted([f for f in os.listdir('.') if f.endswith((".fasta",".fa"))])

    #set subspecies labels to none or parsed from input file
    subspecies = None
    if args.subspecies is not None:
        subspecies = parse_taxa(args.subspecies)
        
    for f in flist:
        if args.relabel == "species":
            label_key = relabel_species(f, fdir, subspecies, args.voucherize)
            write_label_key(label_key, f, ldir)
            
        elif args.relabel == "accession":
            label_key = relabel_accession(f, fdir, subspecies, args.voucherize)
            write_label_key(label_key, f, ldir)
            
        elif args.relabel == "species_acc":
            label_key = relabel_species_acc(f, fdir, subspecies, args.voucherize)
            write_label_key(label_key, f, ldir)
            
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished relabeling alignments. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    

if __name__ == '__main__':
    main()
