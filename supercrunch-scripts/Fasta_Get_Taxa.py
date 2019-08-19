''''
SuperCRUNCH: Fasta_Get_Taxa module

    Fasta_Get_Taxa: Construct 'species' and 'subspecies' label sets from a 
    directory of fasta files containing description lines. Writes two output 
    files which are lists of all unique 'species' and 'subspecies' labels 
    found. For smaller data sets, this can help generate a list of taxa. 
    There are likely to be spurious names in the lists and they should 
    always be examined before using them for other purposes. This is 
    particularly true for the subspecies labels, because records with 
    only binomial names are expected to produce incorrect subspecies labels. 
    There are some filters in place to reduce bad names from being recorded. 
    Although this reduces some of the 'junk' names, there are bound to be 
    some bad names produced.

    In some cases (such as population level data) it may be desirable to
    obtain a name such as 'Genus Species Sample', rather than 'Genus Species
    Subspecies'. In such cases, the 'Sample' component may be a museum code
    or other numerical identifier. This can be obtained by using the optional 
    flag --numerical.
        
    The resulting list files can be quickly edited to produce a combined 
    list of valid species and subspecies names, which can then be used 
    for other steps in SuperCRUNCH.
    
    Input fasta files must have extension '.fasta' or '.fa', to be read.

    Output Files:
    
    Species_Names.txt - List of unique binomial names constructed from 
                        record descriptions. If records are labeled correctly 
                        this should correspond to the genus and species. This 
                        file should be carefully inspected.
                                               
    Subspecies_Names.txt - List of unique trinomial names constructed from 
                           record descriptions. If records actually contain 
                           subspecies labels they will be captured in this list, 
                           however if the records only contain a binomial name 
                           then spurious names may be produced. This file should 
                           be carefully inspected. 
           
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
import argparse
from datetime import datetime

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Fasta_Get_Taxa: Construct 'species' and 'subspecies' label sets from a 
    directory of fasta files containing description lines. Writes two output 
    files which are lists of all unique 'species' and 'subspecies' labels 
    found. For smaller data sets, this can help generate a list of taxa. 
    There are likely to be spurious names in the lists and they should 
    always be examined before using them for other purposes. This is 
    particularly true for the subspecies labels, because records with 
    only binomial names are expected to produce incorrect subspecies labels. 
    There are some filters in place to reduce bad names from being recorded. 
    Although this reduces some of the 'junk' names, there are bound to be 
    some bad names produced.

    In some cases (such as population level data) it may be desirable to
    obtain a name such as 'Genus Species Sample', rather than 'Genus Species
    Subspecies'. In such cases, the 'Sample' component may be a museum code
    or other numerical identifier. This can be obtained by using the optional 
    flag --numerical.
        
    The resulting list files can be quickly edited to produce a combined 
    list of valid species and subspecies names, which can then be used 
    for other steps in SuperCRUNCH.
    
    Input fasta files must have extension '.fasta' or '.fa', to be read.
	---------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory with fasta file(s).")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to write "
                            "output files.")
    
    parser.add_argument("--numerical",
                            action='store_true',
                            help="OPTIONAL: Allows the third part of a name to also contain "
                            "numbers and special characters, rather than just letters. Useful "
                            "for samples of the same species with museum/field codes or other "
                            "unique identifiers immediately following the species label.")
    
    return parser.parse_args()


def avoid_lists():
    """
    This is just a place to store these long lists. They
    can be edited here as necessary. This function returns
    both lists.
    """
    #these are mostly gene related words that generally will lead to spurious names
    genes = ["ALDEHYDE", "AMELOBLASTIN", "AMELOGENIN", "AMELOTIN", "ANDROGEN", "ANTERIOR",
                 "ANTIZYME", "APOLIPOPROTEIN", "ARYL", "BASIC", "CALCIUM", "CONSTITUTIVE",
                 "CYSTEINE", "CYTOCHROME", "DUAL", "ECTODERMAL", "ENAMELIN", "ENDOTHELIN",
                 "EPIDERMAL", "ESTROGEN", "FOLLICLE", "FOLLICULAR", "GLYCOGEN", "GROWTH",
                 "IMMUNOGLOBULIN", "INHIBIN", "INTERFERON", "KERATIN", "KISSPEPTIN",
                 "LEUCINE", "LFNG", "LINKAGE", "LORICRIN", "MEDIATOR", "MEGAKARYOBLASTIC",
                 "MESOTOCIN", "MICRORNA", "MUTS", "NERVE", "NEURAMINIDASE", "NEUROLIGIN",
                 "OPSIN", "ORNITHINE", "OXYTOCIN", "PARAPINOPSIN", "PEPTIDOGLYCAN", "PININ",
                 "PROLACTIN", "PROSTAGLANDIN", "PROTEIN", "PROTOCADHERIN", "PSEUDOGENE",
                 "RECOMBINATION", "RHODOPSIN", "RIBONUCLEASE", "SCAFFOLDIN", "SECRETORY",
                 "SMALL", "SODIUM", "SOLUTE", "SPEXIN", "SYNUCLEIN", "TEASHIRT",
                 "THIOREDOXIN", "THYMIDINE", "TRANSIENT", "TRANSPOSON", "UBINUCLEIN",
                 "VASOTOCIN", "VITELLOGENIN", "ZINC", "INTERNAL", "HISTONE", "NATURAL",
                 "ALPHA", "CARDIAC", "BONE", "OOCYTE", "DYNEIN", "TRANSFER", "GAPDH",
                 "FIBROBLAST", "HEART", "HOMEOBOX", "SONIC", "RING", "FORKHEAD", "PP",
                 "RED", "TWIST", "VERTEBRATE", "ACID", "BETA", "HOMEOBOX", "SACSIN",
                 "STEROIDOGENIC", "TITIN", "SUBUNIT", "PAIRED", "SYNAPSIN", "SYNTAXIN",
                 "PRT", "SRC", "LACTATE", "MYOGLOBIN", "PANCREATIC", "ATPASE", "HEMOGLOBIN",
                 "CONE", "CYCLIC", "DYNEIN", "GUSTDUCIN", "LONG", "PARIETOPSIN", "PHOSPHODIESTERASE",
                 "PINOPSIN", "ROD", "AXONEMAL", "CREATINE", "GENOMIC", "CYTB", "ENOS",
                 "INOS", "NNOS", "CITRATE", "PIN", "PITUITARY", "BRAIN", "GO", "BROTHER",
                 "NEUROTROPHIN", "SHORT", "CSB", "PHOSDUCIN", "OLFACTORY", "TRUNCATED",
                 "LOW", "CMOS", "MYOSIN", "ADRENERGIC", "DOPAMINE", "DOUBLESEX", "GLUTAMATE",
                 "HEAT", "NEURONAL", "PROGESTERONE", "SEROTONIN", "TRYPTOPHAN", "TYROSINE",
                 "PHOSPHOLIPASE", "PUTATIVE", "ASRIN", "COMPLEMENT", "CYSTATIN", "VENOM",
                 "BRADYKININ", "CYSTATIN", "DISINTEGRIN", "HYPOTHETICAL", "METALLOPROTEASE",
                 "PHOSPHOLIPASE", "SERINE", "VASCULAR", "VENOM", "ACETYLCHOLINERG", "DENMOTOXIN",
                 "DOPAMINE", "GLUTAMINYL", "AUSTRALIAN", "GLUTAMINYL", "IRDITOXIN", "POMC",
                 "THREE", "TOXIN", "CANDIDUXIN", "CANDOXIN", "KAPPA", "PHOSPHOLIPAS", "WEAK",
                 "ACETYLCHOLINESTERASE", "ANTIMICROBIAL", "BRAIN", "BUNGARUSKUNIN", "CONTROL",
                 "METALLOPROTEINASE", "MYOSIN", "PHOSPHOLIPASE", "VENOM", "PUTATIVE", "LOCUS",
                 "ETS", "ALDOLASE", "SOLUTION", "EXTRACTION", "AXIS", "BRAIN", "FOS", "KINASE",
                 "PATCHED", "PROSECRETONEURIN", "SECRETED", "TYPE", "THROMBIN", "LYSOSOMAL",
                 "VIMENTIN", "STEM", "FACTOR", "ANONYMOUS", "ARYL", "KINASE", "AMELOGENIN",
                 "BASIC", "BTB", "BASIC", "BONE", "CASPASE", "CARTILAGE", "OOCYTE", "CHEMOKINE",
                 "DYNEIN", "ENDOTHELIN", "ECTODERMAL", "EXOPHILIN", "FOLLICLE", "GALANIN",
                 "GROWTH", "HOLOCARBOXYLASE", "INHIBIN", "SICKLE", "KINESIN", "LEUCINE",
                 "LEUCINE", "MELANOCORTIN", "MEGAKARYOBLASTIC", "LEUKEMIA", "MUTS", "MYOSIN",
                 "NERVE", "NATURAL", "NITRIC", "NEUROTROPHIN", "PHOSDUCIN", "PININ", "PROLACTIN",
                 "PROSTAGLANDIN", "PROTEIN", "RECOMBINATION", "POLYMERASE", "RHODOPSIN", "SOLUTE",
                 "SYNUCLEIN", "SUPRESSOR", "TNF", "UBINUCLEIN", "ZINC", "CYTOCHROME", "DEHYDROGENASE",
                 "UNTRANSLATED", "NEUROTOXIN", "ACIDIC", "PROTHROMBIN", "PROVIRAL", "ULTRA",
                 "HYDROXYSTEROID", "STEROID", "SULFOTRANSFERASE", "THYROID", "TRANSTHYRETIN",
                 "URIDINE", "ARRESTIN", "PHOTORECEPTOR", "VITAMIN", "BLUE", "GREEN", "CHROMODOMAIN",
                 "DOUBLE", "SOLUBLE", "GENOTYPE", "PROGLUCAGON", "BROMODOMAIN", "BUTYRYLCHOLINESTERASE",
                 "CHONDROADHERIN", "DOLICHOL", "GLUCURONIC", "PANNEXIN", "PERIPLAKIN", "PYRUVATE",
                 "SPASTIC", "INTERSPERSED", "FIBRINOGEN", "SERUM", "SATELLITE", "EUKARYOTIC",
                 "SKELETAL", "ATRAGIN", "ATRATOXIN", "BROMODOMAIN", "BUTYRYLCHOLINESTERA",
                 "CARDIOTOXIN", "CHONDROADHERIN", "COBROTOXIN", "DOLICHOL", "GLUCURONIC", "NATRIN",
                 "NATRIURETIC", "PANNEXIN", "PERIPLAKIN", "PYRUVATE", "SPASTIC", "MUSCLE",
                 "CYTOTOXIN", "NEUROTOXIN", "ACIDIC", "NEUTRAL", "UNKNOWN", "ANKYRIN", "CARNITINE",
                 "DELETED", "DERMATAN", "FIBRONECTIN", "FRIZZLED", "FUCOSYLTRANSFERASE", "LIGASE",
                 "MISSING", "POTASSIUM", "RECEPTOR", "SUSHI", "SUPPRESSOR", "STONIN",
                 "PHOSPHOGLUCONATE", "ACETYLCHOLINERGIC", "PREPROGLUCAGON"]

    #these are other words often found in records that will lead to spurious names
    badwords = ["AFF" , "CF" , "SP", "PREDICTED" , "UNVERIFIED" , "UNCULTURED" , "VOUCHER" ,
                    "STRAIN" , "SPECIMEN" , "HAPLOTYPE" , "MITOCHONDRION" , "PARTIAL" , "GENE" ,
                    "MRNA" , "WHOLE" , "CLONE" , "FROM" , "MITOCHONDRIAL" , "ISOLATE" ,
                    "MICROSATELLITE", "CHROMOSOME", "TRANSCRIPTION", "RIBOSOMAL", "X", "COUNTRY",
                    "LIVER", "TRANSCRIPTIONAL", "MITOCHONDRIA", "COMPLETE", "UNPLACED", "TSA"]
        
    return genes, badwords
    
def get_names(f, numerical, genes, badwords):
    """
    Function to extract name labels from fasta file. It will
    examine two-part names (e.g., Genus Species) and three-
    part names (e.g., Genus Species Subspecies/Identifier). 
    There are some basic filters in place to try to reduce 
    the number of spurious names generated.
    """

    #initiate empty sets for two-part (sp) and three-part (ssp) names
    sp_set, ssp_set = set(), set()
    
    #set up accumulator for number of lines read
    rcnt = int(0)
    #open file
    with open(f, 'r') as fh:
        #iterate through lines
        for line in fh:
            #find record labels, denoted by >
            if line.startswith(">"):
                #add to line count and print progress when appropriate
                rcnt += 1
                if rcnt % 10000 == 0:
                    print("\tProcessed {:,} records...".format(rcnt))
                
                #process species (binomial) names here
                #begin by removing punctuation and splitting line
                #and keeping first two elements, if possible
                parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                              if len(line.split()) >= int(3)][1:3]
                
                #if two strings are in the list
                if parts1:
                    #check if string one (Genus) and string two (Species)
                    #are more than a single letter, and only alphabetical
                    if (parts1[0].isalpha()
                            and len(parts1[0]) > 1
                            and parts1[-1].isalpha()
                            and len(parts1[-1]) > 1):
                        #next ensure both strings are not in the badwords list
                        if (parts1[0].upper() not in badwords
                                and parts1[-1].upper() not in badwords):
                            #next ensure both strings are not in the genes list
                            if (parts1[0].upper() not in genes
                                    and parts1[-1].upper() not in genes):
                                #if all checks passed, add to the two-label set
                                taxon_sp = " ".join(parts1)
                                sp_set.add(taxon_sp)

                #process subspecies (trinomial) names here
                #begin by removing punctuation and splitting line
                #and keeping first three elements, if possible
                parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split()
                              if len(line.split()) >= int(4)][1:4]
                #if three strings are in the list
                if parts2:
                    #base checks on whether or not we should allow
                    #the third label to contain numbers and other special characters

                    #here we do not allow any numbers in the third label
                    if not numerical:
                        #check if string one (Genus), string two (Species),
                        #and string three (Subspecies)
                        #are more than a single letter, and only alphabetical
                        #and that string three is not initially all uppercase
                        if (parts2[0].isalpha()
                                and len(parts2[0]) > 1
                                and parts2[1].isalpha()
                                and len(parts2[1]) > 1
                                and parts2[-1].isalpha()
                                and len(parts2[-1]) > 1
                                and not parts2[-1].isupper()):
                            #next ensure all strings are not in the badwords list
                            if (parts2[0].upper() not in badwords
                                    and parts2[1].upper() not in badwords
                                    and parts2[-1].upper() not in badwords):
                                #next ensure all strings are not in the genes list
                                if (parts2[0].upper() not in genes
                                        and parts2[1].upper() not in genes
                                        and parts2[-1].upper() not in genes):
                                    #if all checks passed, add to the three-label set
                                    taxon_ssp = " ".join(parts2)
                                    ssp_set.add(taxon_ssp)

                    #here we do allow numbers in the third label
                    elif numerical:
                        #check if string one (Genus), string two (Species),
                        #are more than a single letter, and only alphabetical
                        #only check if string three (Subspecies/Identifier)
                        #is more than a single character
                        if (parts2[0].isalpha()
                                and len(parts2[0]) > 1
                                and parts2[1].isalpha()
                                and len(parts2[1]) > 1
                                and len(parts2[-1]) > 1):
                            #next ensure all strings are not in the badwords list
                            if (parts2[0].upper() not in badwords
                                    and parts2[1].upper() not in badwords
                                    and parts2[-1].upper() not in badwords):
                                #next ensure all strings are not in the genes list
                                if (parts2[0].upper() not in genes
                                        and parts2[1].upper() not in genes
                                        and parts2[-1].upper() not in genes):
                                    #if all checks passed, add to the three-label set
                                    taxon_ssp = " ".join(parts2)
                                    ssp_set.add(taxon_ssp)
                        
    print("\n\tRead {:,} total records...".format(rcnt))
    						
    return sp_set, ssp_set

def process_names(indir, outdir, numerical):
    """
    Identify fasta files in input directory and obtain
    name labels for each file using the get_names() function.
    The names are added to larger sets and converted to sorted
    lists, which are then written to output files in the output
    directory specified.
    """
    
    os.chdir(indir)
    flist = [f for f in os.listdir('.') if f.endswith(('.fasta', '.fa'))]
    all_sp, all_ssp = set(), set()

    genes, badwords = avoid_lists()
    
    for f in flist:
        print("\nExamining {}:".format(f))
        sp_set, ssp_set = get_names(f, numerical, genes, badwords)
        all_sp.update(sp_set)
        all_ssp.update(ssp_set)
        
    sp, ssp = sorted(all_sp), sorted(all_ssp)

    print("\n\nFound a total of {} unique 'species' names.".format(len(sp)))
    print("Found a total of {} unique 'subspecies' names.\n\n".format(len(ssp)))

    os.chdir(outdir)
    with open("Species_Names.txt", 'a') as fh:
        for s in sp:
            fh.write('{}\n'.format(s))

    with open("Subspecies_Names.txt", 'a') as fh:
        for s in ssp:
            fh.write('{}\n'.format(s))
            
def main():
    tb = datetime.now()
    args = get_args()

    process_names(args.indir, args.outdir, args.numerical)

    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished. Total elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")
    
if __name__ == '__main__':
    main()
