"""
### ======= Align the reads to the composite MeVChiTok genome ======= ###
This snakemake file uses BWA to align the filtered reads to the composite MeVChiTok genome.
These alignments are then sorted and indexed with samtools. Finally, the alignments will be 
used to generate an updated consensus sequence from the brain.  

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""


rule get_gff:
    """ 
    Download the gff format annotation and change the chromosome names.

    This is necessary for our custom, composite reference genome.
    """
    output: join(config['gff_dir'], 'MeVChiTok.gff')
    params: ftp = lambda wildcards: config['MeVChiTok']['gff']
    shell: 
        """
        wget -O - {params.ftp} | gunzip -c | sed 's/NC_001498.1/MeVChiTok/g' > {output}
        """
    

rule bwa_index:
    """ 
    Index the genome with BWA before alignment with BWA. 

    * Note that the reference genome should already be in the directory
    and isn't downloaded from the NCBI FTP site. *
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

        
rule bwa_align:
    """ 
    Perform short read alignment of the filtered 
    viral readas with `bwa-mem`.
    
    Sort the aligned reads with samtools sort.
    """
    input: 
        reads = [join(config['filter_dir'], "{accession}", "{accession}_R1.fastq.gz"), 
                 join(config['filter_dir'], "{accession}", "{accession}_R2.fastq.gz")],
        genome = join(config['index_dir']['bwa'], 'MeVChiTok.fa')
    output: bam = join(config['align_dir'], "{accession}", "{accession}.sorted.bam"),
            bai = join(config['align_dir'], "{accession}", "{accession}.sorted.bam.bai")
    threads: config['threads']['max_cpu']
    params: readgroup = lambda wildcards: get_readgroup(wildcards) 
    conda: '../envs/align.yml'
    shell: 
        """
        bwa mem -t {threads} \
            -R "@RG\\tID:1\\tSM:{params.readgroup}\\tLB:cattaneo\\tPU:illumina" {input.genome} \
            {input.reads} | \
            samtools view -bh | \
            samtools sort -o {output.bam} - 
        samtools index {output.bam}
        """
