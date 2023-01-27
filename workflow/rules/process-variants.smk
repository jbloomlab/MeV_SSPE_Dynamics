"""
### ======= Process variants and assign genotypes ======= ###
This snakemake file runs a series of R notebooks to process the variants and assign haplotypes.

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""

rule quality_control:
    """
    Combine the variants identfied by the variant callers and see how the compare. I 
    also fill in variants that are missed in one tissue by only one caller. 
    """
    input:
        join(config['variant_dir'], "variants.csv")
    output:
        join(config['notebook_dir'], "quality-control.html")
    params: outcsv=join(config['variant_dir'], "filtered_variants.csv")
    conda: "../envs/r.yml"
    script: "../notebooks/quality-control.Rmd"


rule identify_backgrounds:
    """
    I identify the variants that define the genome-1 and genome-2 backgrounds.
    """
    input:
        join(config['notebook_dir'], "quality-control.html")
    output:
        join(config['notebook_dir'], "identify-backgrounds.html")
    params: 
        incsv=join(config['variant_dir'], "filtered_variants.csv"),
        outcsv=join(config['variant_dir'], "genotyped_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
    conda: "../envs/r.yml"
    script: "../notebooks/identify-backgrounds.Rmd"


rule subclonal_haplotypes:
    """
    I identify and number the subclonal haplotypes at high enough frequency with clustering. 
    """
    input:
        join(config['notebook_dir'], "identify-backgrounds.html")
    output:
        join(config['notebook_dir'], "cluster-subclonal-haplotypes.html")
    params: 
        incsv=join(config['variant_dir'], "genotyped_variants.csv"),
        ancestral=join(config['ref_dir'], "annotated_SSPE_consensus_snps.csv"),
        outcsv=join(config['variant_dir'], "clustered_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
    conda: "../envs/r.yml"
    script: "../notebooks/cluster-subclonal-mutations.Rmd"
