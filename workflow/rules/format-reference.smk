### ======= Download and format reference genomes and annotations ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule get_gff:
    """ 
    Download the gff format annotation and change the chromosome names.
    """
    output: join(config['gff_dir'], 'MeVChiTok.gff')
    params: ftp = lambda wildcards: config['MeVChiTok']['gff']
    shell: 
        """
        wget -O - {params.ftp} | gunzip -c > {output}
        sed 's/NC_001498.1/MeVChiTok/g' {output}
        """
    

rule bwa_index:
    """ Index the genome with `BWA` before mapping.
    """
    input: join(config['ref_dir'], 'MeVChiTok.fa')
    output: join(config['index_dir']['bwa'], 'MeVChiTok.fa')
    conda: '../envs/align.yml'
    params: algorithm="bwtsw"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        """


rule samtools_index:
    """ Index genome with `samtools` for `BSQR`. 
    """
    input: join(config['ref_dir'], 'MeVChiTok.fa')
    output: 
        fa=join(config['index_dir']['samtools'], 'MeVChiTok.fa'),
        idx=join(config['index_dir']['samtools'], 'MeVChiTok.fa.fai'),
        idxdict=join(config['index_dir']['samtools'], 'MeVChiTok.fa.dict')
    conda: '../envs/align.yml'
    shell:
        """
        cp {input} {output.fa}
        samtools faidx {output.fa}
        samtools dict -o {output.idxdict} {output.fa}
        """




