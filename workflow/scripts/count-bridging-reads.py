"""
The goal of this notebook is to count the bridging reads between all qualified SNPs.
"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import os
import itertools
from collections import Counter, defaultdict

import pandas as pd 
import numpy as np

import pysam 
from Bio import SeqIO 



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
    

def phase_variants(bampath, snps):
    """
    Get the phase of variants by making a dictionary keyed by 
    read name and add allele counts for all polymorphic sites.
    
    Parameters
    ----------
    
    bampath: str
        Path to a BAM file to phase variants for.
    
    snps: list
        A sorted list of the SNPs to phase.
    
        
    Returns
    -------
    
    pandas.DataFrame
        DataFrame with the pairwise phase of all targeted SNPs
    
    """
    
    # Searchable set of alternative alleles at each position 
    alt_alleles = set((pos, alt) for ref, pos, alt in snps)
    # Searchable set of reference alleles at each position 
    ref_alleles = set((pos, ref) for ref, pos, alt in snps)
    
    # First and last position to visit in the pileup column
    start = snps[0][1]
    stop = snps[-1][1]
    
    # Empty dict to store qnames
    qnames = defaultdict(list)
    # Empty dict to store snp pairs
    phased = defaultdict(list)
    
    ## ==== Go through the bam file and identify haplotypes ==== ##
    
    # Open the alignment file
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        
        # Iterate over the pileup column at each position
        for pileupcolumn in bamfile.pileup(stepper = 'nofilter',
                                           flag_filter = 0,
                                           min_base_quality = 25,
                                           start = start,
                                           stop = stop,):
            
            # Check if the position contains a polymorphic position
            if pileupcolumn.reference_pos + 1 in set(pos for ref, pos, alt in snps):
                
                # Iterate over each read in the pileup column 
                for pileupread in pileupcolumn.pileups:
                    
                    # Check that it's a good read - no deletions, qfails, or skips
                    if check_read(pileupread.alignment) and not pileupread.is_del and not pileupread.is_refskip:

                        # Save the read name
                        qname = pileupread.alignment.query_name

                        # Save the 1-indexed position
                        position = pileupcolumn.pos + 1

                        # Save the base at that position in the read
                        allele = pileupread.alignment.query_sequence[pileupread.query_position]

                        # Check if this read has the SNP or not - assuming bialleleic! 
                        if (position, allele) in alt_alleles:
                            phase = 1 # The read has the alt allele at this position

                        elif (position, allele) in ref_alleles:
                            phase = 0 # The read has the ref allele at this position
                        else:
                            continue 

                        # Add the qname to the dictionary along with phase at position 
                        qnames[qname].append((position, phase, allele))

    ## ==== Collate haplotypes and count for all SNP paris observed ==== ##
    
    for qname, alleles in qnames.items(): 

        # Can't phase anything with a single SNP
        if len(set(alleles)) == 1: 
            continue 

        haplotype = defaultdict(set)
        for pos, phase, allele in alleles:
            haplotype[pos].add(phase)

        # Don't use reads with read pairs that disagree 
        for phases in haplotype.values():
            if len(phases) > 1:
                continue

        # Get a list of allele observations for this read
        allele_obsvs = sorted([(position, list(phase)[0])
                               for position, phase 
                               in haplotype.items()], key = lambda x: x[0])

        # Get the phases of all combinations represented on this read
        for allele_one, allele_two in itertools.combinations(allele_obsvs, 2): 

            allele_pair = (allele_one[0], allele_two[0])

            phasing = f"{allele_one[1]}{allele_two[1]}"

            # Add them to a dictionary indexed by the combination of positions 
            phased[allele_pair].append(phasing)

    # Count the haplotypes for each pair of positions with overlaping reads 
    counts = {comp: Counter(haps) for comp, haps in phased.items()}

    # Convert to a dataframe 
    counts_df = (pd.DataFrame(counts)
                 .T
                 .fillna(0)
                 .reset_index()
                 .rename(columns = {"level_0": "snp_1", "level_1": "snp_2"})
                 .sort_values(by=['snp_1', 'snp_2'])
    )

    # Check that there are no redundant pairs - these would need to be combined 
    assert len(set(
        Counter(
            tuple(sorted((pos_1, pos_2))) 
            for pos_1, pos_2 
            in zip(counts_df.snp_1, counts_df.snp_2))
        .values())) == 1

    return counts_df
                            

def main():

    print("Starting analysis.")

    ## == Paths from Snakemake pipeline == ##

    # Get the path list for the input BAM files from snakemake rule.
    bampath = snakemake.input.bam
    # Dataframe with SNPs and genome-1/2 annotations
    snpspath = str(snakemake.input.snps)
    # Get the path to the output 
    outcsv = str(snakemake.output.csv)


    ## == Main script functionality == ##
    
    # List of all SNPs including the reference allele - assume all sites are biallelic 
    snps_df = pd.read_csv(snpspath)
    snps_list = sorted(
                    list(
                        {(REF, POS, ALT) for REF, POS, ALT in
                        zip(snps_df.REF, snps_df.POS, snps_df.ALT)}
                    ),
                    key = lambda x: x[1])

    
    tissue = " ".join(os.path.basename(bampath).split(".")[0].split("_"))
    print(f'Assigning reads for {tissue}')
    
    phase_df = phase_variants(bampath, snps_list)
    phase_df["Tissue"] = tissue

    phase_df.to_csv(outcsv, index = False)


if __name__ == "__main__":
    main()

