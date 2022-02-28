"""
The goal of this script is to call assign variants in a give tissues to genome-1 or genome-2
"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import os
from collections import Counter, defaultdict
from unicodedata import name
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
    

def assign_reads(genome_dict, bampath):
    """
    Save a list of read names associated with either genome-1 or genome-2 SNPs

    Parameters
    ----------
    genome_dict : dict
        Dictionary of either genome-1 or genome-2 only SNPs.
        
    bampath: str
        Path to bamfile to assign reads for. 

    Returns
    -------
    dict
        Dictionary of read name sets assigned to either genome-1 or genome-2 key. 
    """
    
    # Assigning reads to either genome 1 or genome 2
    assert len(genome_dict.keys()) == 2, "There must be two genotypes."
    assert "genome-1" in genome_dict.keys() and "genome-2" in genome_dict.keys(), "Must be genome-1 and genome-2."
            
    # Empty list to store qnames
    qnames = defaultdict(list)
    
    # Open the alignment file
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        # Iterate through the bams and get reads with genome 1 or genome 2 SNPs 
        for genotype, snps in genome_dict.items():
            # Iterate over the pileup column at each position
            for pileupcolumn in bamfile.pileup(stepper = 'nofilter', flag_filter = 0, min_base_quality = 25):
                # Check if the position contains a genotypeable position
                if pileupcolumn.reference_pos + 1 in set(pos for pos, snp in snps):
                    for pileupread in pileupcolumn.pileups:
                        # Check that it's a good read
                        if check_read(pileupread.alignment) and not pileupread.is_del and not pileupread.is_refskip:
                            # Save the qname
                            if (pileupcolumn.reference_pos + 1, pileupread.alignment.query_sequence[pileupread.query_position]) in snps:
                                qnames[genotype].append(pileupread.alignment.query_name)
                            else:
                                if genotype == "genome-1":
                                    qnames["genome-2"].append(pileupread.alignment.query_name)
                                else:
                                    qnames["genome-1"].append(pileupread.alignment.query_name) 
                                
    return {k:set(v) for k,v in qnames.items()}


def subset_bam(bampath, outpath, qnames):
    """
    Take a bamfile and subset it into a new bam file 
    based on a list of read names.
    
    Parameters
    ----------
    bampath : str
        Path pointing to an indexed Bam format file.
        
    read_names: list
        A list of readnames to extract the Bam file. 
        
    prefix: str
        The prefix of the file to be written out during the analysis.
        
    """
    with pysam.AlignmentFile(bampath) as bam:
        with pysam.AlignmentFile(outpath, "w", template = bam) as outbam:
            for read in bam.fetch(until_eof=True):
                if read.query_name in qnames:
                    outbam.write(read)
     # Message about status
    print(f"Wrote bam subset to {outpath}\n\nDone!\n")
    

def sort_bam(bampath):
    """
    Take a bamfile and sort and index it.
    
    Parameters
    ----------
    bampath : str
        Path pointing to an indexed Bam format file.

    """
    
    # Sort the bam 
    print("Sorting Bam\n")
    pysam.sort("-o", f"{os.path.splitext(bampath)[0]}.sorted.bam", bampath)
    
    # Index the sorted bam
    print("Indexing Bam\n")
    pysam.index(f"{os.path.splitext(bampath)[0]}.sorted.bam")   
    
    print("Done!\n")
        

def count_coverage(bampath,
                   snps,
                   callback_function = check_read, 
                   contig = "MeVChiTok",
                   ref_path = "../../config/ref/MeVChiTok-SSPE.fa", 
                   minimum_qual = 25):
    """
    Count the coverage of all SNPs of interest in a filtered BAM file. 
    
    Parameters
    ----------
    bampath : str
        Path pointing to an indexed Bam format file.

    snps : list
        List of all SNPs to filter resulting dataframe by
        
    Returns
    -------
    pd.DataFrame
        Dataframe of SNP counts and total coverage in an alignment file. 
        
    """
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        
        # The count_coverage method counts the occurances of each base at each position. 
        # It excludes reads based on the callback function.
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT",
                                                                               bamfile.count_coverage(contig=contig,
                                                                               read_callback=callback_function,
                                                                               quality_threshold=minimum_qual))})
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        # Add the position (1-indexed)
        count_df['POS'] = count_df.index + 1
        # Add the reference allele
        count_df['REF'] = [base.upper() for base in list(SeqIO.parse(ref_path, "fasta"))[0].seq]
        # convert counts to frequency 
        count_df.iloc[:,0:4] = count_df.iloc[:,0:4].div(count_df.DP, axis = 0)
        # handle any NaNs created by dividing by 0 coverage
        count_df = count_df.fillna(0)
        # Melt the data frame to a longer ('tidy') form
        count_df = pd.melt(count_df, 
                           id_vars=['POS', 'DP', 'REF'],
                           value_vars=[base for base in 'ATGC'],
                           value_name='AF',
                           var_name='ALT')
        # Remove anything with 0 coverage.
        count_df = count_df[count_df['AF'] > 0]
        # TRUE/FALSE if it's a SNP
        count_df['IS_SNP'] = np.where(count_df['ALT'] != count_df['REF'], True, False)
        # Add the number of times a given allele is observed.
        count_df['OBSV'] = count_df.DP * count_df.AF
        
        # Filter to only get SNPs
        count_df.loc[count_df['IS_SNP']]
        count_df = count_df.drop(columns=['IS_SNP'])
        
        # New colunm of the SNP to filter by
        count_df['SNP'] = count_df["REF"] + count_df["POS"].astype(str) + count_df["ALT"]
        
        return count_df[count_df['SNP'].isin(snps)]
        

def main():

    print("Starting analysis.")

    ## == Paths from Snakemake pipeline == ##
    # Get the path list for the input BAM files from snakemake rule.
    bampath = snakemake.input.bam
    # Dataframe with SNPs and genome-1/2 annotations
    snpspath = str(snakemake.input.snps)
    # Get the path to the reference genome.
    refpath = str(snakemake.input.genome)
    # Get the path to the output 
    outdir = str(snakemake.output.dir)
    # Get the path to the final dataframe
    outcsv = str(snakemake.output.csv)
    # Get the contig
    contig = str(snakemake.params.contig)

    ## == Main script functionality == ##
    
    # Import the dataframe with mutations labled by identity
    genomes_df = pd.read_csv(snpspath)

    print("Sucessfully imported labeled genomes.")

    # Isolate the genome-1 and genome-2 only mutations. 
    genome_dict = {
        'genome-1': set(zip(genomes_df.loc[genomes_df.Haplotype == "genome-1"].POS, genomes_df.loc[genomes_df.Haplotype == "genome-1"].ALT)),
        'genome-2': set(zip(genomes_df.loc[genomes_df.Haplotype == "genome-2"].POS, genomes_df.loc[genomes_df.Haplotype == "genome-2"].ALT))
    }

    # Get the name of this tissue to add to the dataframe.
    tissue = " ".join(snakemake.wildcards.accession.split("_"))

    print(f"Starting read splitting for {tissue}")


     # Get the reads from a bam belonging to each tissue. 
    qnames_dict = assign_reads(genome_dict, bampath)

    # How many reads get assigned to both genomes (could be recombinats, reversions, or template switching)
    reads_in_both = qnames_dict['genome-1'] & qnames_dict['genome-2']

    print(f"{len(reads_in_both)/len(qnames_dict['genome-1'] | qnames_dict['genome-2']) * 100} percent of reads are in both sets of genomes.\n")

    genotype_df_list = []
    # For each genotype, get the count of SNPs
    for genotype, qnames in qnames_dict.items():
        qnames = qnames - reads_in_both
        # Subset the bam for only these reads
        subset_bam(bampath, os.path.join(outdir, f"{genotype}.bam"), qnames)
        # Sort and index for variant id
        sort_bam(os.path.join(outdir, f"{genotype}.bam"))
        # Count the coverage for SNPs of interest 
        coverage_df = count_coverage(os.path.join(outdir, f"{genotype}.sorted.bam"),
                                     snps=genomes_df.SNP.to_list(),
                                     callback_function = check_read, 
                                     contig = contig,
                                     ref_path = refpath, 
                                     minimum_qual = 25)

        coverage_df['genotype'] = f"{genotype}"
        genotype_df_list.append(coverage_df)

    genotype_df = pd.concat(genotype_df_list)
    genotype_df['Tissue'] = f"{tissue}"

    print(f"Writing dataframe for {tissue}")

    # Write out the final dataframe
    genotype_df.to_csv(outcsv, index = False)
    
if __name__ == "__main__":
    main()

