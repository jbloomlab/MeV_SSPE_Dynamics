### ======= Preprocess the fastq files ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule fetch_fastq:
    """ 
    Move the fastq files into the results directory and compress them. 
    """
    output: R1 = join(config['fastq_dir'], "{accession}", "{accession}_R1.fastq.gz"),
            R2 = join(config['fastq_dir'], "{accession}", "{accession}_R2.fastq.gz")
    params: R1_path = lambda wildcards: SAMPLE_DF.loc[(SAMPLE_DF.Run == wildcards.accession)].R1.item(),
            R2_path = lambda wildcards: SAMPLE_DF.loc[(SAMPLE_DF.Run == wildcards.accession)].R2.item()
    shell: "gzip -c {params.R1_path} > {output.R1} && gzip -c {params.R2_path} > {output.R2}"


rule trim_adapters:
    """
    Fast all-in-one processing of single `fastq` files: 
    1. Automatic adaptor trimming
    2. Low-qual base filtering (50% of bases w/ Phred >20) 
    3. Reporting by HTML and JSON.
    """
    input:  
        R1 = join(config['fastq_dir'], "{accession}", "{accession}_R1.fastq.gz"),
        R2 = join(config['fastq_dir'], "{accession}", "{accession}_R2.fastq.gz")
    output:
        R1 = join(config['trim_dir'], "{accession}", "{accession}_R1.fastq.gz"),
        R2 = join(config['trim_dir'], "{accession}", "{accession}_R2.fastq.gz"),
        html = join(config['trim_dir'], "{accession}", "{accession}.fastp.html"),
        json = join(config['trim_dir'], "{accession}", "{accession}.fastp.json")
    conda: '../envs/preprocessing.yml'
    shell:
        """ 
        fastp \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
            --html {output.html} \
            --json {output.json} 
        """


rule filter_reads:
    """ 
    Filter out reads using kmer filtering w/ BBduk Paried-End implementation.
    This will only work for the viral genomes because larger genomes
    begin to require huge amounts of heap memory.

    Genome is determined from the `samples.csv` file configuration.
    The genomes are indexed with samtools. 
    """
    input: 
        R1 = join(config['trim_dir'], "{accession}", "{accession}_R1.fastq.gz"),
        R2 = join(config['trim_dir'], "{accession}", "{accession}_R2.fastq.gz"),
        genome = join(config['ref_dir'], "MeVChiTok.fa")
    output: 
        R1 = join(config['filter_dir'], "{accession}", "{accession}_R1.fastq.gz"), 
        R2 = join(config['filter_dir'], "{accession}", "{accession}_R2.fastq.gz"),
        R1_unmatched = join(config['filter_dir'], "{accession}", "{accession}_R1.unmatched.fastq.gz"), 
        R2_unmatched = join(config['filter_dir'], "{accession}", "{accession}_R2.unmatched.fastq.gz"),
        stats = join(config['filter_dir'], "{accession}", "{accession}.filter.stats")
    threads: config['threads']['max_cpu']
    params: error = join(config['filter_dir'], "{accession}", "{accession}.error.log")
    conda: '../envs/preprocessing.yml'
    shell:
        """
        bbduk.sh -Xmx80g \
            in1={input.R1} \
            in2={input.R2} \
            out1={output.R1_unmatched} \
            out2={output.R2_unmatched} \
            outm1={output.R1} \
            outm2={output.R2} \
            ref={input.genome} \
            k=31 \
            hdist=2 \
            stats={output.stats} \
            overwrite=TRUE \
            t={threads} \
            &> {params.error}
        """
