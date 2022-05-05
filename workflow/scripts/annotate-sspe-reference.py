"""
The goal of this script is to annotate the coding effect of SSPE consensus mutations.

The SSPE consensus only includes the mutations that are fixed in *most* of the SSPE samples. 
"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2022 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import os 
import pandas as pd
import numpy as np
from Bio import SeqIO

def parse_gff(gff_path, genes_to_keep = ['N', 'P', 'M', 'F', 'H', 'L']):
    """
    Parse the Measles GFF file to get the coordinates of the coding sequences.
    """
    # Hold the parsed genes
    gff_dict = {}
    
    # Parse the GFF file
    with open(gff_path) as gff_file:
        for line in gff_file:
            # Ignore the header
            if line.startswith('#'):
                continue
            # Only take the coding sequences
            else:
                record = line.strip().split("\t")
                if record[2] == "CDS":
                   
                    # Get the start and stop 
                    start = int(record[3])
                    stop = int(record[4])
                    
                    # Get the gene and product names
                    annot_dict = {annot.split("=")[0]: annot.split("=")[1] for annot in record[8].split(";")}
                    gene = annot_dict['gene']

                    # There are four P/V/C reading frames
                    if gene == 'P/V/C':
                        gene = annot_dict['product'][0].upper()
                    
                    # Only keep the annotations of interest
                    if gene in genes_to_keep:
                        gff_dict[gene] = [start, stop]
                        
    return gff_dict


def parse_reference(cds_positions, reference):
    """
    Parse the reference sequence into annotated coding sequence dictionary.
    """
    # Save the annotated coding sequences in a dictionary
    cds_dict = {}
    
    # Get the codon sequence and position for every gene
    for gene, index in cds_positions.items():

        cds = reference[(index[0]-1):(index[1])]

        cds_dict[gene] = [(i, (i + index[0]), nt) for i, nt in enumerate(cds)]
        
    return cds_dict


def get_codon_index(index):
    
    if index % 3 == 0:
        return 1
    elif index % 3 == 1:
        return 2
    elif index % 3 == 2:
        return 3


def translate(codon):
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    codon : str
        three letter DNA sequence

    Returns
    -------
    str
        one letter amino acid code

    Raises
    ------
    AssertionError
        error if codon sequence is invalid
        
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    
    codon = "".join(codon)
    
    assert codon in table.keys(), "Not a valid codon sequence."
    
    return table[codon]

def main():

    # Path to original reference
    refpath = str(snakemake.input.genome)
    # Path to measles GFF
    gffpath = str(snakemake.input.gff)
    # Path to SSPE consensus mutations
    csvpath = str(snakemake.input.csv)
    
    # Path to output annotated consensus mutations.
    outcsv = str(snakemake.output)

    # Consensus mutations as dataframe
    consensus_snps_df = pd.read_csv(csvpath)
    # Reference sequence as list
    reference_seq = [base for record in SeqIO.parse(refpath, "fasta") for base in record.seq]
    # Dictionary of conding regions
    cds_dict = parse_reference(parse_gff(gffpath), reference_seq)
    # List of the position and alt base for SSPE consensus SNPs
    consensus_snps = [(POS, ALT) for POS, ALT in zip(consensus_snps_df.POS, consensus_snps_df.ALT)]

    # Determine the coding effect of each conensus mutation
    annotated_snps = list()
    for pos, alt in consensus_snps:
        
        for gene, indicies in parse_gff(gffpath).items():
            
            if pos >= indicies[0] and pos <= indicies[1]:
                gene_name = gene
                gene_pos = cds_dict[gene][pos-indicies[0]][0]
                codon_pos = get_codon_index(gene_pos)
                break
            else:
                gene_name = "intergenic"
                
        if gene_name != "intergenic":

            if codon_pos == 1:
                codon = [cds_dict[gene][pos-indicies[0]][2],
                        cds_dict[gene][pos+1-indicies[0]][2],
                        cds_dict[gene][pos+2-indicies[0]][2]]

                codon_alt = codon.copy()
                codon_alt[codon_pos-1] = alt

                coding_pos = int((gene_pos+3)/3)
                ref = codon_alt[codon_pos-1]

            elif codon_pos == 2:
                codon = [cds_dict[gene][pos-1-indicies[0]][2],
                        cds_dict[gene][pos-indicies[0]][2],
                        cds_dict[gene][pos+1-indicies[0]][2]]

                codon_alt = codon.copy()
                codon_alt[codon_pos-1] = alt

                coding_pos = int((gene_pos+2)/3)
                ref = codon_alt[codon_pos-1]

            elif codon_pos == 3:
                codon = [cds_dict[gene][pos-2-indicies[0]][2],
                        cds_dict[gene][pos-1-indicies[0]][2],
                        cds_dict[gene][pos-indicies[0]][2]]

                codon_alt = codon.copy()
                codon_alt[codon_pos-1] = alt

                coding_pos = int((gene_pos+1)/3)
                ref = codon_alt[codon_pos-1]

            annotated_snps.append((pos, reference_seq[pos-1], alt, gene_name, translate(codon), translate(codon_alt), coding_pos))
            
        else: 
        
            annotated_snps.append((pos, reference_seq[pos-1], alt, gene_name, np.nan, np.nan, np.nan))
                
    annotated_df = pd.DataFrame(list(set(annotated_snps)), columns = ["POS", "REF", "ALT", "Gene", "WT_AA", "MUT_AA", "POS_AA"])
    annotated_df.to_csv(outcsv, index=False)


if __name__ == '__main__':
    main()


