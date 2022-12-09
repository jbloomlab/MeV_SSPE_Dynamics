"""
### ======= Make SSPE reference and realign filtered reads ======= ###
This snakemake file makes an SSPE reference from the consensus mutations 
identified in SSPE samples aligned against the combined MeVChiTok reference.

After making the SSPE reference, the filtered reads are realigned to this reference.

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""

rule make_sspe_reference: 
    """
    Make a reference that incorporates fixed mutations.
    
    We're considering a mutation to be fixed in the brain if
    it's in *most* of the SSPE samples, currently 12/15 at more
    than 90% frequency in each tissue sample.
    """
    input: bam = expand(join(config['align_dir'], "{accession}", "{accession}.sorted.bam"), accession = samples),
           genome = join(config['ref_dir'], 'MeVChiTok.fa')
    output: fasta = join(config['ref_dir'], 'MeVChiTok-SSPE.fa'),
            csv = join(config['ref_dir'], 'SSPE_consensus_snps.csv')
    params: contig = "MeVChiTok",
            nobsv = 12, 
            minfreq = .90
    threads: config['threads']['max_cpu']
    conda: '../envs/pysam.yml'
    script: '../scripts/make-sspe-reference.py'


rule bwa_sspe_index:
    """ 
    Index the SSPE genome with BWA before realignment with BWA. 
    """
    input: join(config['ref_dir'], 'MeVChiTok-SSPE.fa')
    output: join(config['index_dir']['bwa'], 'MeVChiTok-SSPE.fa')
    conda: '../envs/align.yml'
    params: algorithm="bwtsw"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        """


rule samtools_sspe_index:
    """ 
    Index the SSPE genome with `samtools` for `BSQR`. 

    This index is important for variant calling and other
    downstream analyses that require random access to the genome.
    """
    input: join(config['ref_dir'], 'MeVChiTok-SSPE.fa')
    output: 
        fa=join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa'),
        idx=join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa.fai'),
        idxdict=join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa.dict')
    conda: '../envs/align.yml'
    shell:
        """
        cp {input} {output.fa}
        samtools faidx {output.fa}
        samtools dict -o {output.idxdict} {output.fa}
        """


rule annotate_consensus_mutations:
    """
    Since we're curious about the coding effect of some of these mutations,
    this script annotates the mutations with the coding effect.
    """
    input: 
        csv=join(config['ref_dir'], 'SSPE_consensus_snps.csv'),
        gff=join(config['gff_dir'], 'MeVChiTok.gff'),
        genome=join(config['ref_dir'], 'MeVChiTok.fa')
    output: join(config['ref_dir'], 'annotated_SSPE_consensus_snps.csv')
    conda: '../envs/pysam.yml'
    script: '../scripts/annotate-sspe-reference.py'


rule bwa_realign:
    """ 
    Perform short read alignment of the filtered viral reads with `bwa-mem`.
    This time, realigning to SSPE consensus generated above.
    
    Sort the aligned reads with samtools sort and index the sorted bam file.
    """
    input: 
        reads = [join(config['filter_dir'], "{accession}", "{accession}_R1.fastq.gz"), 
                 join(config['filter_dir'], "{accession}", "{accession}_R2.fastq.gz")],
        genome = join(config['index_dir']['bwa'], 'MeVChiTok-SSPE.fa')
    output: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam"),
            bai = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam.bai")
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
