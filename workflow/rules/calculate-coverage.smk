### ======= Calculate read coverage ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule samtools_depth:
    """ Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam")
    output: join(config['coverage_dir'], "{accession}", "{accession}.depth")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -m 0 -a -q {params.score} -g DUP {input.bam} > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{accession}", "{accession}.depth"), accession=samples)
    output: depth = join(config['coverage_dir'], "merged.depth"),
            header = temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "POS\tDP\tAccession"}}1' {output.header} > {output.depth}
        """