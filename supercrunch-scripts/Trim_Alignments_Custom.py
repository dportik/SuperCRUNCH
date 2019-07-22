''''
SuperCRUNCH: Trim_Alignments_Custom module

	Trim_Alignments_Custom: Use a customizable trimming function to trim all
    alignments. There are three main options (-a), which include edge trimming 
    and internal trimming (edges_internal), edge trimming only (edges), and 
    internal trimming only (internal). The trimming for each of these options
    is guided by the optional arguments (-w, -p, -t, -d, and -l). These arguments
    have default values which can be changed by adding the relevant flags with
    new values. After custom trimming, the alignments are 'cleaned up' using 
    trimAl.
    
    The output format must be specified using the -f argument (choices: fasta, nexus, 
    or phylip). Timmed alignment output files will be moved to following output 
    directory specified. 

    In general, this trimming function is very aggressive in removing columns, 
    particularly when sequence length heterogeneity is relatively high.
    An output file called 'Trimming_Summary.txt' is created, which has information 
    about the starting lengths and trimmed lengths for every input alignment. If
    too many columns are being removed, you can adjust the optional arguments and
    inspect this output file to see the effect on the number of columns removed, 
    and adjust accordingly.

    Several functions present in this customizable trimming module are 
    modified versions of functions included in the PHYLUCE package, specifically 
    the generic_align.py module.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

-------------------------
Compatible with Python 2.7 & 3.7
Dependencies:
	-trimal (installed in path)
-------------------------

SuperCRUNCH project
https://github.com/dportik/SuperCRUNCH
Written by Daniel Portik 
daniel.portik@gmail.com
July 2019
Distributed under the 
GNU General Public Lincense
'''

import os
import re
import numpy
import argparse
import shutil
import subprocess as sp
from datetime import datetime
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""-----------------------------------------------------------------------------
	Trim_Alignments_Custom: Use a customizable trimming function to trim all
    alignments. There are three main options (-a), which include edge trimming 
    and internal trimming (edges_internal), edge trimming only (edges), and 
    internal trimming only (internal). The trimming for each of these options
    is guided by the optional arguments (-w, -p, -t, -d, and -l). These arguments
    have default values which can be changed by adding the relevant flags with
    new values. After custom trimming, the alignments are 'cleaned up' using 
    trimAl.
    
    The output format must be specified using the -f argument (choices: fasta, nexus, 
    or phylip). Timmed alignment output files will be moved to following output 
    directory specified. 

    In general, this trimming function is very aggressive in removing columns, 
    particularly when sequence length heterogeneity is relatively high.
    An output file called 'Trimming_Summary.txt' is created, which has information 
    about the starting lengths and trimmed lengths for every input alignment. If
    too many columns are being removed, you can adjust the optional arguments and
    inspect this output file to see the effect on the number of columns removed, 
    and adjust accordingly.

    Several functions present in this customizable trimming module are 
    modified versions of functions included in the PHYLUCE package, specifically 
    the generic_align.py module.

    Input fasta files should be labeled as 'NAME.fasta' or 'NAME.fa', 
    where NAME represents the gene/locus. The NAME portion should not 
    contain any periods or spaces, but can contain underscores. Output 
    files are labeled using a prefix identical to NAME.

    DEPENDENCIES: Executables in path: trimal.
	-----------------------------------------------------------------------------""")
    
    parser.add_argument("-i", "--indir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains "
                            "the aligned and relabeled fasta files to trim.")
    
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-a", "--analysis",
                            required=True,
                            choices=["edges_internal", "edges", "internal"],
                            help="REQUIRED: Specify the method for trimming alignments.")
    
    parser.add_argument("-f", "--format",
                            required=True,
                            choices=["fasta", "nexus", "phylip"],
                            help="REQUIRED: Specify the output file format for trimmed "
                            "alignment.")
    
    parser.add_argument("-w", "--window",
                            required=False,
                            type=int,
                            default=20,
                            help="Optional: Sliding window size for trimming.")
    
    parser.add_argument("-p", "--proportion",
                            required=False,
                            type=float,
                            default=0.65,
                            help="Optional: The proportion of taxa required to have "
                            "sequence at alignment ends.")
    
    parser.add_argument("-t", "--threshold",
                            required=False,
                            type=float,
                            default=0.65,
                            help="Optional: The proportion of residues required across "
                            "the window in proportion of taxa.")
    
    parser.add_argument("-d", "--max_divergence",
                            required=False,
                            type=float,
                            default=0.20,
                            help="Optional: The max proportion of sequence divergence allowed "
                            "between any row of the alignment and the consensus sequence.")
    
    parser.add_argument("-l", "--min_length",
                            required=False,
                            type=int,
                            default=100,
                            help="Optional: The minimum length of alignments to keep.")
    
    return parser.parse_args()


def record_formatter(trim, name):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    return a string formatted as a biopython sequence record
    """
    return SeqRecord(Seq(trim, Gapped(IUPAC.ambiguous_dna, "-?")), id=name, name=name, description=name)
    
def alignment_consensus(alignment):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    return consensus for an alignment object using BioPython
    """
    consensus = []
    for pos in xrange(alignment.get_alignment_length()):
        col = alignment[:,pos].upper()
        cnt = Counter(col)
        base, occurrence = cnt.most_common(1)[0]
        proportion = float(occurrence)/len(col)
        if proportion >= 0.5:
            consensus.append(base)
        else:
            consensus.append('N')
            
    return ''.join(consensus)
    
def get_ends(seq):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    Find the start and end of sequence data for a given alignment row
    """
    f = re.compile("^([-]+)")
    result = f.search(str(seq.seq))
    if result:
        start_gap = len(result.groups()[0])
    else:
        start_gap = 0
    r = re.compile("([-]+)$")
    result = r.search(str(seq.seq))
    if result:
        end_gap = len(result.groups()[0])
    else:
        end_gap = 0
        
    return start_gap, len(seq.seq) - end_gap

def running_average(alignment, window_size, proportion, threshold):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    Trim an alignment (assuming default parameters) such that `proportion`
    of taxa have `threshold` of residues at each position across the first
    and last `window_size` column slice of the alignment
    """
    good_alignment = []
    taxa = len(alignment)
    majority_of_characters = int(round(proportion * taxa, 0))
    for column in xrange(alignment.get_alignment_length()):
        column_count = Counter(alignment[:, column])
        if column_count['-'] <= majority_of_characters:
            del column_count['-']
            if column_count.most_common(1)[0][1] >= majority_of_characters:
                good_alignment.append(True)
            else:
                good_alignment.append(False)
        else:
            good_alignment.append(False)
    good_alignment = numpy.array(good_alignment)
    for start_clip in xrange(good_alignment.size):
        if good_alignment[start_clip] != False:
            window = good_alignment[start_clip: start_clip + window_size]
            proportion = float(sum(window)) / len(window)
            if proportion > threshold:
                break
    reverse_good_alignment = good_alignment[::-1]
    for end_clip in xrange(reverse_good_alignment.size):
        if reverse_good_alignment[end_clip] != False:
            window = reverse_good_alignment[end_clip: end_clip + window_size]
            proportion = float(sum(window)) / len(window)
            if proportion >= threshold:
                end_clip = reverse_good_alignment.size - end_clip
                break
            
    return start_clip, end_clip

def stage_one_trimming(alignment, window_size, proportion, threshold, min_len):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    First stage alignment trimming to find and trim edges of a given
    alignment.  Calls running_average function above to determine reasonable
    alignment start and end trimming for the entire alignment block.
    """
    start, end = running_average(alignment, window_size, proportion, threshold)
    s1_trimmed = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna, "-?"))
    for sequence in alignment:
        sequence.seq.alphabet = IUPAC.IUPACAmbiguousDNA()
        if start >= 0 and end:
            trim = sequence[start:end]
            if set(trim) != set(['-']) and set(trim) != (['?']) and len(trim) >= min_len:
                s1_trimmed.append(sequence[start:end])
            else:
                s1_trimmed = None
                break
        else:
            s1_trimmed = None
            break
        
    return s1_trimmed

def stage_two_trimming(s1_trimmed, window_size, max_divergence, min_len):
    """
    ---------------------------------------------------------------------
    MODIFIED FUNCTION FROM PHYLUCE: generic_align.py
    ---------------------------------------------------------------------
    Alignment row-by-row trimming.  After stage one trimming, iterate
    over rows of alignment to find differences between the alignment
    consensus and the row (taxon) of data.  Trim those ends that differ
    from the consensus with > `divergence` across a `window_size` window.
    Goes to third round of filtering to remove edges that end up with only '----'
    characters to start or end alignment block.
    """
    s2_trimmed = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna, "-?"))
    consensus_array = numpy.array(list(alignment_consensus(s1_trimmed)))
    for sequence in s1_trimmed:
        sequence = sequence.upper()
        start, end = get_ends(sequence)
        orig_seq_array = numpy.array(list(sequence))
        seq_array = orig_seq_array[start:end]
        bad_start = 0
        bad_end = len(sequence)
        compare = (seq_array != consensus_array[start:end])
        for bad_start in xrange(compare.size):
            window = compare[bad_start: bad_start + window_size]
            divergence = float(sum(window))/window.size
            if divergence < max_divergence:
                break
        reversed_compare = compare[::-1]
        for bad_end in xrange(reversed_compare.size):
            window = reversed_compare[bad_end: bad_end + window_size]
            divergence = float(sum(window))/window.size
            if divergence < max_divergence:
                bad_end = reversed_compare.size - bad_end
                break
        orig_seq_array[:start + bad_start] = '-'
        orig_seq_array[start + bad_end:] = '-'
        trim = ''.join(orig_seq_array)
        if set(trim) != set(['-']) and set(trim) != (['?']) and len(trim) >= min_len:
            s2_trimmed.append(record_formatter(trim, sequence.id))
        else:
            s2_trimmed = None
            break
        
    return s2_trimmed

def trim_alignment(aln, analysis, window_size=20, proportion=0.65,
                       threshold=0.65, max_divergence=0.20, min_len=100):
    """
    Trim a given alignment using 1) trim alignment block ends and trim row data, 
    2) trim alignment block ends only, or 3) trim row data only. Drops alignments
    shorter than 100 bp.
    """
    
    if analysis == "edges_internal":
        s1_trimmed = stage_one_trimming(aln, window_size, proportion, threshold, min_len)
        
        if s1_trimmed:
            trimmed = stage_two_trimming(s1_trimmed, window_size, max_divergence, min_len)
        else:
            trimmed = None
            
    elif analysis == "edges":
        trimmed = stage_one_trimming(aln, window_size, proportion, threshold, min_len)
        
    elif analysis == "internal":
        trimmed = stage_two_trimming(aln, window_size, max_divergence, min_len)
        
    return trimmed

def run_trim(f, analysis, outformat, logname, window, proportion,
                 threshold, max_divergence, min_length):
    """
    Run trimming method on alignment. Afterwards, eliminate bad columns 
    using gap-threshold value of 0.05 and any columns composed
    of only gaps using trimal -noallgaps, and output in format desired. 
    Write a log file of the starting length and final length of the
    alignments.
    """
    #read in alignment using biopython
    alignment = AlignIO.read(f, 'fasta')
    print("\nTrimming {}".format(f))
    
    #get length of starting alignment
    slen = alignment.get_alignment_length()

    #run trimming routine using modified phyluce functions
    trimmed = trim_alignment(alignment, analysis, window,
                                 proportion, threshold, max_divergence, min_length)
    
    #if alignment was successfully trimmed, do things
    if trimmed:
        #write trimmed alignment to temporary file
        #so that we can do final cleanup with trimal
        prefix = f.split('.')[0]
        temp = "{}.temp.fasta".format(prefix)
        with open(temp, 'w') as fh:
            fh.write(trimmed.format('fasta'))
        
        #command for trimal gap-threshold trimming with value of 0.05
        #output in fasta because will transform this file again
        tempout = "{}_trimmed.temp.fasta".format(prefix)
        cmd1 = ("trimal -in {0} -out {1} -fasta -gt 0.05"
                    .format(temp, tempout))
        #execute command
        proc = sp.call(cmd1, shell=True)

        #generate command for trimal noallgaps feature, but output
        #in format selected by user - this is the final trimmed file
        if outformat == "fasta":
            outname = "{}_trimmed.fasta".format(prefix)
            cmd2 = ("trimal -in {0} -out {1} -fasta -noallgaps"
                        .format(tempout, outname))
            
        elif outformat == "nexus":
            outname = "{}_trimmed.nex".format(prefix)
            cmd2 = ("trimal -in {0} -out {1} -nexus -noallgaps"
                        .format(tempout, outname))
                
        elif outformat == "phylip":
            outname = "{}_trimmed.phy".format(prefix)
            cmd2 = ("trimal -in {0} -out {1} -phylip_paml -noallgaps"
                        .format(tempout, outname))
        #execute command
        proc = sp.call(cmd2, shell=True)

        #read final trimmed alignment using biopython
        final_aln = AlignIO.read(outname, outformat)
        #obtain the length of the final alignment
        flen = final_aln.get_alignment_length()
        
        #print information to screen and write to log file
        print("\tLength of starting alignment: {}.".format(slen))
        print("\tLength of trimmed alignment: {}.".format(flen))
        with open(logname, 'a') as fh:
            fh.write("{}\t{}\t{}\t{}\n".format(f, slen, flen, (int(slen)-int(flen))))
        
    #if alignment was overly trimmed and failed, it will
    #return None - handle this case here
    else:
        print("\n\t***Alignment did not meet requirements after trimming.***\n")
        with open(logname, 'a') as fh:
            fh.write("{}\t{}\t{}\t{}\n".format(f, slen, "**FAILED**", "**NA**"))

def make_dirs(outdir, outformat):
    """
    Creates directory path names and makes output directories.
    Returns directory paths, which are used to move around  
    output files during cleanup steps.
    """
    os.chdir(outdir)
    #get current path
    curpath = os.getcwd()
            
    #create paths using os.path.join() to avoid any issues
    tdir = os.path.join(curpath, "Trimmed-{}-Files".format(outformat))
    if not os.path.exists(tdir):
        os.mkdir(tdir)

    sdir = os.path.join(curpath, "Summary-File".format(outformat))
    if not os.path.exists(sdir):
        os.mkdir(sdir)

    return tdir, sdir
    

def cleanup(outdir):
    """
    Moves relevant output files to their output directories. 
    Deletes temporary or empty files.
    """
    trimmed = [shutil.move(f, outdir) for f in os.listdir('.')
                   if f.endswith(("_trimmed.fasta", "_trimmed.nex", "_trimmed.phy"))]
    trimmed = [os.remove(f) for f in os.listdir('.') if f.endswith(".temp.fasta")]
    
def write_log(log_list, outdir):
    """
    Simple function to take contents of log_list
    and write to output file.
    """
    logname = "Trimming_Summary.txt"
    
    with open(logname, 'a') as fh:
        fh.write("Fasta\tStart_Length\tTrimmed_Length\tBases_Trimmed\n")
        
    with open(logname, 'a') as fh:
        for l in log_list:
            fh.write("{}\t{}\t{}\t{}\n".format(l[0], l[1], l[2], (int(l[1])-int(l[2]))))
    
def main():
    args = get_args()
    tb = datetime.now()
    tdir, sdir = make_dirs(args.outdir, args.format)
    
    os.chdir(args.indir)
    fastas = sorted([f for f in os.listdir('.') if f.endswith((".fa", ".fasta"))])
    
    logname = "Trimming_Summary.txt"
    with open(logname, 'a') as fh:
        fh.write("Fasta\tStart_Length\tTrimmed_Length\tBases_Trimmed\n")    

    for f in fastas:
        b = datetime.now()
        run_trim(f, args.analysis, args.format, logname, args.window,
                     args.proportion, args.threshold, args.max_divergence, args.min_length)
        cleanup(tdir)
        f = datetime.now()
        e = f - b
        print("\tElapsed time: {0} (H:M:S)\n".format(e))
    shutil.move(logname, sdir)
    
    tf = datetime.now()
    te = tf - tb
    print("\n\n--------------------------------------------------------------------------------------")
    print("\nFinished trimming alignments. Elapsed time: {0} (H:M:S)\n".format(te))
    print("--------------------------------------------------------------------------------------\n\n")    
    
if __name__ == '__main__':
    main()
