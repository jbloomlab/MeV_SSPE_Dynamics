### ======= Identify haplotypes with Clique-SNV ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule cliqueSNV:
    input: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam")
    output: sam = temp(join(config['haplotype_dir'], "{accession}", "{accession}.sorted.sam")),
            outdir = directory(join(config['haplotype_dir'], "{accession}"))
    params:
        prefix = join(config['haplotype_dir'], "{accession}", "{accession}.sorted"),
        cliqueSNV = join(config['tools'], 'clique-snv.v2.0.3.jar')
    conda: "../envs/haplotype.yml"
    threads: config['threads']['max_cpu']
    shell: 
        """
        # Convert the BAM file to SAM
        samtools view -h -o {output.sam} {input.bam}

        # Run cliqueSNV on the SAM w/ VCF output
        java -Xmx80g -jar {params.cliqueSNV} \
            -m snv-illumina \
            -in {output.sam} \
            -ignoreDeletion \
            -t 40 \
            -tf 0.02 \
            -el 700 \
            -sp 56 \
            -ep 15786 \
            -os 56 \
            -oe 15786 \
            -threads {threads} \
            -outDir {output.outdir} \
            -log 
        """