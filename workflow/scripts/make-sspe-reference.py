"""
The goal of this script is to make an SSPE consensus sequence for a variant calling reference. 

The SSPE consensus will only include the mutations that are fixed in most of the SSPE samples. 
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
import numpy as np
import os 
from collections import Counter


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
                   minimum_QUAL = 25): 
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

    
def call_variants(bampath,
                  refpath,
                  contig,
                  threads,
                  minimum_QUAL = 25,
                  minimum_COV = 100):
    """
    Read in BAM file and convert to a dataframe containing the frequency of 
    of any bases present at a given position in the reference genome using the
    `pysam` command `count_coverage`. 

    Parameters
    ----------
    filepath : str
        path to the bam file to be parsed
        
    callback_function : function
        function that decides which reads to keep/exclude
        
    ref : str
        name of the contig to count coverage over
        
    ref_path : str
        path to the reference genome as fasta
         
    minimum_qual : int
        minimum QUAL score to count base observation at a position.

    Returns
    -------
    Pandas.DataFrame
       Data Frame containing the bases represented at each positon in the genome
        
    """
    
    # Regions for parallelization
    regions = split(get_ref_sequence(refpath), threads)

    # Multithreadding to make the count dataframe 
    with concurrent.futures.ProcessPoolExecutor() as executor: 

        results = executor.map(functools.partial(count_coverage,
                                                 bampath,
                                                 contig,
                                                 refpath,  
                                                 minimum_QUAL = minimum_QUAL), 
                               regions)

        count_df = pd.concat([result for result in results]).reset_index()
    
    

    # Melt the data frame to a longer ('tidy') form
    count_df = pd.melt(count_df, 
                       id_vars=['POS', 'DP', 'REF'],
                       value_vars=[base for base in 'ATGC'],
                       value_name='AF',
                       var_name='ALT')
    
    # Remove anything with 0 coverage.
    count_df = count_df[count_df['AF'] > 0]
    # TRUE/FALSE if it's a SNP
    count_df['SNP'] = np.where(count_df['ALT'] != count_df['REF'], True, False)
    # TRUE/FALSE if it's a consensus base.
    count_df['CONS'] = count_df['AF'].map(lambda x: x >= 0.5)
    # Add the number of times a given allele is observed.
    count_df['OBSV'] = count_df.DP * count_df.AF
    # Filter to only get SNPs
    count_df = count_df.loc[count_df['SNP']]

    # Merge the pileup dict and the snp dict. 
    return count_df


def update_reference(bampaths, 
                     refpath,
                     outpath,
                     contig, 
                     threads,
                     fixed_threshold = .95,
                     n_obsv = 11,
                     minimum_QUAL = 25,
                     minimum_COV = 100):
    """
    Update the current MevChiTok reference sequence by incorporating fixed mutations present
    in the alignments to the original reference.

    Parameters
    ----------
    bampaths : list
        a list of paths to the bam file to be parsed
        
    refpath : str
        path to the reference genome as fasta     
       
    outpath : path
        path to the updated reference fasta
        
    contig : str
        name of the contig to count coverage over
        
    threads : int
        number of threads to call variants in parallel
        
    fixed_threshold : int
        minimum AF to be considered a fixed mutation.

    n_obsv: int 
        minimum number of samples to be considered fixed in whole brain.
         
    minimum_QUAL : int
        minimum QUAL score to count base observation at a position.
        
    minimum_COV : int
        minimum DP score to count base observation at a position.

    """
    
    # Go through and call variants in each bam file.
    variant_df_list = []
    for bampath in bampaths: 
        # Call the variants from the alignment
        variant_df = call_variants(bampath,
                              refpath,
                              contig,
                              threads,
                              minimum_QUAL = 25,
                              minimum_COV = 100)
        # Recover the fixed variants 
        variant_df = variant_df[variant_df['AF'] >= fixed_threshold]
        # Append the sample name
        variant_df["Tissue"] = " ".join(os.path.basename(bampath).split(".")[0].split("_"))
        # Append the variant df to the list
        variant_df_list.append(variant_df)
    
    # Make a single dataframe with every fixed variant
    all_fixed_variants_df = pd.concat(variant_df_list)
    
    # Get the fixed SNPs and their total count in the samples 
    fixed_snps = Counter((pos, alt) for pos, alt in zip(all_fixed_variants_df.POS, all_fixed_variants_df.ALT))
    
    # Get the SNPs deemed to be in enough samples to be included in the updated reference
    consensus_snps = [snp for snp, count in fixed_snps.items() if count >= n_obsv]
    
    # Update the reference sequence 
    ref_seq = get_ref_sequence(refpath)

    for pos, base in consensus_snps: 
        ref_seq[pos-1] = base

    # Write the updated reference out to a fasta file. 
    with open(outpath, "w") as outfile:
        outfile.write(">" + "MeVChiTok" + "\n" + "".join(ref_seq) + "\n")
        
    print("Done!")


def main():

    # Get the path list for the input BAM files from snakemake rule.
    bampaths = snakemake.input.bam
    # Get the path to the reference genome.
    refpath = str(snakemake.input.genome)
    # Get the path to the output multi-fasta file. 
    outpath = str(snakemake.output)
    # Get the contig
    contig = str(snakemake.params.contig)
    # Get the number of threads to use
    threads = int(snakemake.threads)

    update_reference(bampaths, 
                    refpath,
                    outpath,
                    contig, 
                    threads,
                    fixed_threshold = .95,
                    n_obsv = 12,
                    minimum_QUAL = 25,
                    minimum_COV = 100)

if __name__ == '__main__':
    main()


