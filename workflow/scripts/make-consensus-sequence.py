"""
The goal of this script is to make a consensus sequence for each tissue. 

The output of this script is used to call variants relative to the consensus.
"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

from Bio import SeqIO #reading fasta format
import pysam #count variant alleles from BAM
import pandas as pd #data frames
import concurrent.futures #multithreadding manager
import functools #handle function arguments in map


def get_ref_sequence(path):
    """
    Function imports the reference sequence as a list of bases. 

    Parameters
    ----------
    path : str
        Path to reference sequence file - should be a fasta file. 

    Returns
    -------
    list
        List of bases in the reference sequence.
        
    """
    return [base.upper() for base in list(SeqIO.parse(path, "fasta"))[0].seq]


def split(refseq, threads):
    """
    Function split reference into approxamatley equal sized chunks
    based on the number of threads alloted.

    Parameters
    ----------
    refseq : list
       Reference sequence.
       
    threads: int
        Integer to the number of threads alloted. 

    Returns
    -------
    list
        List of tuples containing the start and end position for coverage. 
        
    """
    
    k, m = divmod(len(refseq), threads)
    
    return [( i * k + min(i, m), (i+1) * k + min(i+1, m)) for i in range(threads)]


def check_read(read):
    """
    Helper function to decide what reads should
    be keep when parsing alignment file with `pysam`. 

    Parameters
    ----------
    read : AlignedSegment
        read from alignment file parsed with `pysam`.

    Returns
    -------
    bool
        True/False if read should be included
        
    """
    # Exclude Quality Failures
    if read.is_qcfail:
        return False
    # Exclude Secondary Mappings
    if read.is_secondary:
        return False
    # Exclude Unmapped Reads
    if read.is_unmapped:
        return False
    else:
        return True
    

def count_coverage(bampath,
                   contig, 
                   refpath, 
                   region,
                   callback_function = check_read, 
                   minimum_QUAL = 25, 
                   minimum_COV = 100): 
    """
    Helper function to decide what reads should
    be keep when parsing alignment file with `pysam`. 

    Parameters
    ----------
    read : AlignedSegment
        read from alignment file parsed with `pysam`.

    Returns
    -------
    bool
        True/False if read should be included
        
    """

    # Open alignment with pysam
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        # Get a dataframe of the counts
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT",
                                                                               bamfile.count_coverage(contig = contig,
                                                                                                      start = region[0], 
                                                                                                      stop = region[1],
                                                                                                      read_callback = callback_function,
                                                                                                      quality_threshold = minimum_QUAL))})
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        # Add the position 
        count_df['POS'] = [pos + 1 for pos in range(region[0], region[1])]
        # Add the reference allele
        count_df['REF'] = get_ref_sequence(refpath)[region[0]:region[1]]
        # Convert counts to frequency 
        count_df.iloc[:,0:4] = count_df.iloc[:,0:4].div(count_df.DP, axis = 0)
        # Handle any NaNs created by dividing by 0 coverage
        count_df = count_df.fillna(0)
    
    return count_df

    
def make_consensus(bampath,
                   refpath,
                   contig,
                   outpath,
                   threads,
                   minimum_QUAL = 25,
                   minimum_COV = 100):
    """
    Take the count_df, which is output from `count_coverage`, and make a consensus seq as a string. 
    
    Parameters for consensus calling are the minimum quality score, coverage, and the reads that
    can be included (marked duplicates, quality fails, etc..) defined by `check_read`. 
    
    Returns consensus as string
    """
        
    # String to hold the final consensus seq
    consensus = ""
    
    # Regions for parallelization
    regions = split(get_ref_sequence(refpath), threads)

    # Multithreadding to make the count dataframe 
    with concurrent.futures.ProcessPoolExecutor() as executor: 

        results = executor.map(functools.partial(count_coverage,
                                                 bampath,
                                                 contig,
                                                 refpath,  
                                                 minimum_QUAL = minimum_QUAL, 
                                                 minimum_COV = minimum_COV), 
                               regions)

        count_df = pd.concat([result for result in results]).reset_index()
    
    # Iterate over each row and determine the conesnus sequence and build consensus string
    for index, row in count_df.iterrows(): 
        if row.DP < minimum_COV:
            consensus += row.REF
        else: 
            consensus += max([(row["A"], "A"), (row["C"], "C"), (row["G"], "G"), (row["T"], "T")],key=lambda x:x[0])[1]

    
    # Write out a fasta file for phylogenetic analysis
    with open(outpath, "w") as outfile:
        outfile.write(">" + contig + "\n" + consensus + "\n")


def main(): 
    """
    Main function to run the analysis.
    """

    # Get the path list for the input BAM files from snakemake rule.
    bampath = snakemake.input.bam
    # Get the path to the reference genome.
    refpath = str(snakemake.input.genome)
    # Get the path to the output multi-fasta file. 
    outpath = str(snakemake.output)
    # Get the contig
    contig = str(snakemake.params.contig)
    # Get the number of threads to use
    threads = int(snakemake.threads)
    
    # Make the consensus sequence 
    make_consensus(bampath = bampath,
                   refpath = refpath,
                   contig = contig,
                   outpath = outpath,
                   threads = threads,
                   minimum_QUAL = 25,
                   minimum_COV = 100)


if __name__ == '__main__': 
    main()

