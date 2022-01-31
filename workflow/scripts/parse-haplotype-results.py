"""
The goal of this parse the haplotype results from CliqueSNV into a dataframe. 

"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import json
import pandas as pd
from Bio import SeqIO

# Inputs from snakemake 
jsonpath = str(snakemake.input.json)

fastapath = str(snakemake.input.fasta)

refpath = str(snakemake.input.genome)

outpath = str(snakemake.output)

# Function to get reference sequence  
get_ref = lambda path: [base.upper() for base in list(SeqIO.parse(path, "fasta"))[0].seq]

# Read the parameters from json file into a dictionary  
with open(jsonpath) as jsonfile:
    data = json.load(jsonfile)

# Read Fasta of haplotypes into a dictionary
haplotypes = {record.id : record.seq for record in SeqIO.parse(fastapath, "fasta")}

# Parameters for the haplotyping 
t = float(data['settings']['-t'])
tf = float(data['settings']['-tf'])
start = int(data['settings']['-sp'])
stop = int(data['settings']['-ep'])
tissue = " ".join(snakemake.wildcards.accession.split("_"))

# Get SNPs from the haplotypes relative to the reference 
rows = []
for header, haplotype in haplotypes.items():
    # Get the haplotype number and frequency from the record ID
    haplotype_number = header.split("_")[0]
    haplotype_freq = round(float(header.split("_")[-1]), ndigits = 4)
    
    for i, bases in enumerate(zip(get_ref(refpath)[start:stop], [base for base in haplotype])):
        # Get the SNPs relative to the reference
        if len(set(bases)) != 1: 
            rows.append([(i + 1 + start), bases[0], bases[1], haplotype_freq, haplotype_number, t, tf, start, stop, tissue])
            
# Make the dataframe 
outdf = pd.DataFrame(rows, columns=["POS", "REF", "ALT", "HF", "HAPID", "T", "TF", "Start", "Stop", "Tissue"])
outdf.to_csv(outpath, index = False)


