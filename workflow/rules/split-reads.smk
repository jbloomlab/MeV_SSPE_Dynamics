### ======= Spatial Analysis and Read Splitting ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule genotype_snps:
    """
    Genotype SNPs by counting there occurence on either genome-1 or genome-2 reads.
    """
    input: 
        bam=join(config['realign_dir'], "{accession}", "{accession}.sorted.bam"),
        snps=join(config['spatial_dir'], "labeled_variants.csv"),
        genome=join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa')
    output: 
        dir=directory(join(config['split_dir'], "{accession}")),
        csv=join(config['split_dir'], "{accession}", "{accession}.genotyped.csv")
    params:
        contig="MeVChiTok"
    conda: "../envs/pysam.yml"
    script: "../scripts/split-reads-by-genotype.py"

rule aggregate_genotype:
    """
    This rule aggregates all of the genotype counts. 
    """
    input: expand(join(config['split_dir'], "{accession}",  "{accession}.genotyped.csv"), accession = samples)
    output: join(config['split_dir'], "genotyped.csv")
    run: aggregate_csv(input, output) 