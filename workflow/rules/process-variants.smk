### ======= Process variants and assign genotypes ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule quality_control:
    input:
        join(config['variant_dir'], "variants.csv")
    output:
        join(config['notebook_dir'], "quality-control.html")
    params: outcsv=join(config['variant_dir'], "filtered_variants.csv")
    conda: "../envs/r.yml"
    script: "../notebooks/quality-control.Rmd"


    