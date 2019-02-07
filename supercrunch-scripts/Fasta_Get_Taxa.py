''''
SuperCRUNCH: Fasta_Get_Taxa module

Usage: python Fasta_Convert.py -i [directory with all fasta files]

    Fasta_Get_Taxa: Construct 'species' and 'subspecies' label sets from a directory of 
    fasta files containing description lines. Writes two output files which are lists of
    all unique 'species' and 'subspecies' labels found. For smaller data sets, this can 
    help generate a list of taxa. There are likely to be spurious names in the
    lists and they should always be examined before using them for other purposes. This is
    particularly true for the subspecies labels, because records with only binomial names
    are expected to produce incorrect subspecies labels. There are some filters in place to
    reduce bad names from being recorded. Although this reduces some of the 'junk' names, there 
    are bound to be other bad names produced.
        
    The resulting list files can be quickly edited to produce a combined list of valid 
    species and subspecies names, which can then be used for other steps in SuperCRUNCH.
    
    Empirical fasta files must have extension '.fasta' or '.fa', to be read.

    Output Files:
    
    Species_Names.txt - List of unique binomial names constructed from record descriptions.
                        If records are labeled correctly this should correspond to the 
                        genus and species. This file should be carefully inspected.
                                               
    Subspecies_Names.txt - List of unique trinomial names constructed from record descriptions.
                           If records actually contain subspecies labels they will be captured
                           in this list, however if the records only contain a binomial name then 
                           spurious names may be produced. This file should be carefully inspected. 
           
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
import argparse
from datetime import datetime

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""---------------------------------------------------------------------------
    Fasta_Get_Taxa: Construct 'species' and 'subspecies' label sets from a directory of 
    fasta files containing description lines. Writes two output files which are lists of
    all unique 'species' and 'subspecies' labels found. For smaller data sets, this can 
    help generate a list of taxa. There are likely to be spurious names in the
    lists and they should always be examined before using them for other purposes. This is
    particularly true for the subspecies labels, because records with only binomial names
    are expected to produce incorrect subspecies labels. There are some filters in place to
    reduce bad names from being recorded. Although this reduces some of the 'junk' names, there 
    are bound to be other bad names produced.
    
    The resulting list files can be quickly edited to produce a combined list of valid 
    species and subspecies names, which can then be used for other steps in SuperCRUNCH.
    
    Empirical fasta files must have extension '.fasta' or '.fa', to be read.

    Output Files:
    
    Species_Names.txt - List of unique binomial names constructed from record descriptions.
                        If records are labeled correctly this should correspond to the 
                        genus and species. This file should be carefully inspected.
                                               
    Subspecies_Names.txt - List of unique trinomial names constructed from record descriptions.
                           If records actually contain subspecies labels they will be captured
                           in this list, however if the records only contain a binomial name then 
                           spurious names may be produced. This file should be carefully inspected. 
	---------------------------------------------------------------------------""")
    parser.add_argument("-i", "--input", required=True, help="REQUIRED: The full path to a directory of fasta files.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: The full path to an existing directory to write output files.")
    return parser.parse_args()

def get_taxa(f):
    '''
    Function to get accession numbers.
    '''
    genes = ["ALDEHYDE", "AMELOBLASTIN", "AMELOGENIN", "AMELOTIN", "ANDROGEN", "ANTERIOR", "ANTIZYME", "APOLIPOPROTEIN", "ARYL", "BASIC", "CALCIUM", "CONSTITUTIVE", "CYSTEINE", "CYTOCHROME", "DUAL", "ECTODERMAL", "ENAMELIN", "ENDOTHELIN", "EPIDERMAL", "ESTROGEN", "FOLLICLE", "FOLLICULAR", "GLYCOGEN", "GROWTH", "IMMUNOGLOBULIN", "INHIBIN", "INTERFERON", "KERATIN", "KISSPEPTIN", "LEUCINE", "LFNG", "LINKAGE", "LORICRIN", "MEDIATOR", "MEGAKARYOBLASTIC", "MESOTOCIN", "MICRORNA", "MUTS", "NERVE", "NEURAMINIDASE", "NEUROLIGIN", "OPSIN", "ORNITHINE", "OXYTOCIN", "PARAPINOPSIN", "PEPTIDOGLYCAN", "PININ", "PROLACTIN", "PROSTAGLANDIN", "PROTEIN", "PROTOCADHERIN", "PSEUDOGENE", "RECOMBINATION", "RHODOPSIN", "RIBONUCLEASE", "SCAFFOLDIN", "SECRETORY", "SMALL", "SODIUM", "SOLUTE", "SPEXIN", "SYNUCLEIN", "TEASHIRT", "THIOREDOXIN", "THYMIDINE", "TRANSIENT", "TRANSPOSON", "UBINUCLEIN", "VASOTOCIN", "VITELLOGENIN", "ZINC", "INTERNAL", "HISTONE", "NATURAL", "ALPHA", "CARDIAC", "BONE", "OOCYTE", "DYNEIN", "TRANSFER", "GAPDH", "FIBROBLAST", "HEART", "HOMEOBOX", "SONIC", "RING", "FORKHEAD", "PP", "RED", "TWIST", "VERTEBRATE", "ACID", "BETA", "HOMEOBOX", "SACSIN", "STEROIDOGENIC", "TITIN", "SUBUNIT", "PAIRED", "SYNAPSIN", "SYNTAXIN", "PRT", "SRC", "LACTATE", "MYOGLOBIN", "PANCREATIC", "ATPASE", "HEMOGLOBIN", "CONE", "CYCLIC", "DYNEIN", "GUSTDUCIN", "LONG", "PARIETOPSIN", "PHOSPHODIESTERASE", "PINOPSIN", "ROD", "AXONEMAL", "CREATINE", "GENOMIC", "CYTB", "ENOS", "INOS", "NNOS", "CITRATE", "PIN", "PITUITARY", "BRAIN", "GO", "BROTHER", "NEUROTROPHIN", "SHORT", "CSB", "PHOSDUCIN", "OLFACTORY", "TRUNCATED", "LOW", "CMOS", "MYOSIN", "ADRENERGIC", "DOPAMINE", "DOUBLESEX", "GLUTAMATE", "HEAT", "NEURONAL", "PROGESTERONE", "SEROTONIN", "TRYPTOPHAN", "TYROSINE", "PHOSPHOLIPASE", "PUTATIVE", "ASRIN", "COMPLEMENT", "CYSTATIN", "VENOM", "BRADYKININ", "CYSTATIN", "DISINTEGRIN", "HYPOTHETICAL", "METALLOPROTEASE", "PHOSPHOLIPASE", "SERINE", "VASCULAR", "VENOM", "ACETYLCHOLINERG", "DENMOTOXIN", "DOPAMINE", "GLUTAMINYL", "AUSTRALIAN", "GLUTAMINYL", "IRDITOXIN", "POMC", "THREE", "TOXIN", "CANDIDUXIN", "CANDOXIN", "KAPPA", "PHOSPHOLIPAS", "WEAK", "ACETYLCHOLINESTERASE", "ANTIMICROBIAL", "BRAIN", "BUNGARUSKUNIN", "CONTROL", "METALLOPROTEINASE", "MYOSIN", "PHOSPHOLIPASE", "VENOM", "PUTATIVE", "LOCUS", "ETS", "ALDOLASE", "SOLUTION", "EXTRACTION", "AXIS", "BRAIN", "FOS", "KINASE", "PATCHED", "PROSECRETONEURIN", "SECRETED", "TYPE", "THROMBIN", "LYSOSOMAL", "VIMENTIN", "STEM", "FACTOR", "ANONYMOUS", "ARYL", "KINASE", "AMELOGENIN", "BASIC", "BTB", "BASIC", "BONE", "CASPASE", "CARTILAGE", "OOCYTE", "CHEMOKINE", "DYNEIN", "ENDOTHELIN", "ECTODERMAL", "EXOPHILIN", "FOLLICLE", "GALANIN", "GROWTH", "HOLOCARBOXYLASE", "INHIBIN", "SICKLE", "KINESIN", "LEUCINE", "LEUCINE", "MELANOCORTIN", "MEGAKARYOBLASTIC", "LEUKEMIA", "MUTS", "MYOSIN", "NERVE", "NATURAL", "NITRIC", "NEUROTROPHIN", "PHOSDUCIN", "PININ", "PROLACTIN", "PROSTAGLANDIN", "PROTEIN", "RECOMBINATION", "POLYMERASE", "RHODOPSIN", "SOLUTE", "SYNUCLEIN", "SUPRESSOR", "TNF", "UBINUCLEIN", "ZINC", "CYTOCHROME", "DEHYDROGENASE", "UNTRANSLATED", "NEUROTOXIN", "ACIDIC", "PROTHROMBIN", "PROVIRAL", "ULTRA", "HYDROXYSTEROID", "STEROID", "SULFOTRANSFERASE", "THYROID", "TRANSTHYRETIN", "URIDINE", "ARRESTIN", "PHOTORECEPTOR", "VITAMIN", "BLUE", "GREEN", "CHROMODOMAIN", "DOUBLE", "SOLUBLE", "GENOTYPE", "PROGLUCAGON", "BROMODOMAIN", "BUTYRYLCHOLINESTERASE", "CHONDROADHERIN", "DOLICHOL", "GLUCURONIC", "PANNEXIN", "PERIPLAKIN", "PYRUVATE", "SPASTIC", "INTERSPERSED", "FIBRINOGEN", "SERUM", "SATELLITE", "EUKARYOTIC", "SKELETAL", "ATRAGIN", "ATRATOXIN", "BROMODOMAIN", "BUTYRYLCHOLINESTERA", "CARDIOTOXIN", "CHONDROADHERIN", "COBROTOXIN", "DOLICHOL", "GLUCURONIC", "NATRIN", "NATRIURETIC", "PANNEXIN", "PERIPLAKIN", "PYRUVATE", "SPASTIC", "MUSCLE", "CYTOTOXIN", "NEUROTOXIN", "ACIDIC", "NEUTRAL", "UNKNOWN", "ANKYRIN", "CARNITINE", "DELETED", "DERMATAN", "FIBRONECTIN", "FRIZZLED", "FUCOSYLTRANSFERASE", "LIGASE", "MISSING", "POTASSIUM", "RECEPTOR", "SUSHI", "SUPPRESSOR", "STONIN", "PHOSPHOGLUCONATE", "ACETYLCHOLINERGIC", "PREPROGLUCAGON"]
    badwords = ["AFF" , "CF" , "SP", "PREDICTED" , "UNVERIFIED" , "UNCULTURED" , "VOUCHER" , "STRAIN" , "SPECIMEN" , "HAPLOTYPE" , "MITOCHONDRION" , "PARTIAL" , "GENE" , "MRNA" , "WHOLE" , "CLONE" , "FROM" , "MITOCHONDRIAL" , "ISOLATE" , "MICROSATELLITE", "CHROMOSOME", "TRANSCRIPTION", "RIBOSOMAL", "X", "COUNTRY", "LIVER", "TRANSCRIPTIONAL", "MITOCHONDRIA", "COMPLETE", "UNPLACED", "TSA"]

    sp_set = set()
    ssp_set = set()

    rcnt = int(0)
    with open(f, 'r') as fh:
    	for line in fh:
    		rcnt += 1
    		if rcnt % 10000 == 0:
    			print "\tRead {} lines...".format(rcnt)
    			
        	if line.startswith(">"):
        		parts1 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if len(line.split()) >= int(3)][1:3]
        		if parts1:
        			if parts1[0].isalpha() and len(parts1[0]) > 1 and parts1[-1].isalpha() and len(parts1[-1]) > 1:
						if parts1[0].upper() not in badwords and parts1[-1].upper() not in badwords:
							if parts1[0].upper() not in genes and parts1[-1].upper() not in genes:
								taxon_sp = " ".join(parts1)
								sp_set.add(taxon_sp)
    						
    			parts2 = [l.replace(",",'').replace(";",'').replace(":",'') for l in line.split() if len(line.split()) >= int(4)][1:4]
    			if parts2:
    				if parts2[0].isalpha() and len(parts2[0]) > 1 and parts2[1].isalpha() and len(parts2[1]) > 1 and parts2[-1].isalpha() and len(parts2[-1]) > 1 and not parts2[-1].isupper():
						if parts2[0].upper() not in badwords and parts2[1].upper() not in badwords and parts2[-1].upper() not in badwords:
							if parts2[0].upper() not in genes and parts2[1].upper() not in genes and parts2[-1].upper() not in genes:
								taxon_ssp = " ".join(parts2)
								ssp_set.add(taxon_ssp)
								
	print "\tRead {} total lines...".format(rcnt)
    						
    return sp_set, ssp_set

def write_accs(dir, out_dir):
    os.chdir(dir)
    f_list = [f for f in os.listdir('.') if f.endswith('.fa') or f.endswith('.fasta')]
    all_sp = set()
    all_ssp = set()
    for f in f_list:
    	print "\nExamining {}\n".format(f)
        sp_set, ssp_set = get_taxa(f)
        all_sp.update(sp_set)
        all_ssp.update(ssp_set)
    sp = sorted(all_sp)
    ssp = sorted(all_ssp)
    print "\n\nFound a total of {} unique 'species' names.".format(len(sp))
    print "Found a total of {} unique 'subspecies' names.\n\n".format(len(ssp))
    
    os.chdir(out_dir)
    with open("Species_Names.txt", 'a') as fh_out:
        for s in sp:
            fh_out.write('{}\n'.format(s))
    with open("Subspecies_Names.txt", 'a') as fh_out:
        for s in ssp:
            fh_out.write('{}\n'.format(s))
def main():
	tb = datetime.now()
	args = get_args()
	write_accs(args.input, args.out_dir)
	tf = datetime.now()
	te = tf - tb
	print "\n\nFinished. Total elapsed time: {0} (H:M:S)\n\n".format(te)
    
    
if __name__ == '__main__':
    main()
