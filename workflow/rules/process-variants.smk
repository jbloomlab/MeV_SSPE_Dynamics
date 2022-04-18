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


rule identify_backgrounds:
    input:
        join(config['variant_dir'], "filtered_variants.csv"),
        join(config['notebook_dir'], "quality-control.html")
    output:
        join(config['notebook_dir'], "identify-backgrounds.html")
    params: 
        outcsv=join(config['variant_dir'], "genotyped_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
        samples=config['samples']['file']
    conda: "../envs/r.yml"
    script: "../notebooks/identify-backgrounds.Rmd"

