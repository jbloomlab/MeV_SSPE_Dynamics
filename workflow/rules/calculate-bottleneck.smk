### ======= Calculate the bottleneck between people ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule calculate_tissue_bottleneck: 
    """
    This rule calculates the 'transmission' bottleneck between tissuses.
    """
    input: join(config['spatial_dir'], "labeled_subhaplotypes.csv")
    output: join(config['spatial_dir'], "{accession}",  "{accession}.bottleneck.csv")
    conda: "../envs/r.yml"
    params: max_bottleneck = 200
    script: "../scripts/calculate-tissue-bottleneck.R"


rule aggregate_bottleneck:
    """
    This rule aggregates all of the bottlenecks. 
    """
    input: expand(join(config['spatial_dir'], "{accession}",  "{accession}.bottleneck.csv"), accession = samples)
    output: join(config['spatial_dir'], "bottleneck.csv")
    run: aggregate_csv(input, output) 
