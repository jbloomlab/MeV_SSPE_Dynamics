"""
### ======= Pilot Analysis of SSPE 1 & 2 brain samples ======= ###
This snakemake file uses the MeVChiTok genome as a reference and calls variants 
on the SSPE 1 and SSPE 2 brain samples before being realigned to the patient specific 
reference. 

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""

rule samtools_index:
    """ 
    Index the SSPE genome with `samtools` for `BSQR`. 

    This index is important for variant calling and other
    downstream analyses that require random access to the genome.
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


rule pilot_samtools_mpileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.

    This is necessary for variant calling with VarScan.
    """
    input: bam = join(config['align_dir'], "{sample}", "{sample}.sorted.bam"),
           genome = join(config['index_dir']['samtools'], 'MeVChiTok.fa')
    output: join(config['pilot_dir'], "{sample}.mpileup.txt")
    params: score = config['BQ'],
            maxdepth = config['maxdepth']
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d {params.maxdepth} -E --excl-flags UNMAP,SECONDARY,QCFAIL -q {params.score} -Q {params.score} -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"


rule pilot_varscan_calling:
    """ 
    SNP calling with Varscan.
    
    Parameters are controlled from the config file.  
    """
    input: 
        varscan = join(config['tools'], "VarScan.v2.4.0.jar"),
        pileup = join(config['pilot_dir'], "{sample}.mpileup.txt")
    output: 
        variants = join(config['pilot_dir'], "{sample}.varscan.vcf")
    params:
        minimum_coverage = config['min_coverage'],
        minumum_supporting_reads = config['min_reads_supporting'],
        minimum_base_quality = config['BQ'],
        minimum_variant_freq = config['min_allele_frequency'],
        strand_filter = config['strand_bias_filter']
    conda: '../envs/java.yml'    
    shell:
        """
        # Call SNPs and InDels using the mpileup file
        java -jar {input.varscan} \
            mpileup2cns {input.pileup} \
            --variants 1 \
            --output-vcf 1 \
            --min-coverage {params.minimum_coverage} \
            --min-reads2 {params.minumum_supporting_reads} \
            --min-avg-qual {params.minimum_base_quality} \
            --strand-filter {params.strand_filter} \
            --min-var-freq {params.minimum_variant_freq} > {output.variants}
        """

rule pilot_vcf_to_table:
    """
    Convert the VCF files to tables for easy data
    analysis in R or Python.
    """
    input: join(config['pilot_dir'], "{sample}.varscan.vcf")
    output: join(config['pilot_dir'], "{sample}.varscan.txt")
    conda: '../envs/gatk.yml'    
    shell:
        """
        if [ -s {input} ]; then
            gatk VariantsToTable \
                -V {input} \
                -F CHROM -F POS -F QUAL -F REF -F ALT \
                -F DP -F AF -F FILTER -GF DP \
                -GF RD -GF FREQ -GF SDP -GF AD -F ANN \
                -O {output}
        else
            touch {output}
        fi
        """


rule pilot_add_metadata:
    """
    This rule adds in the metadata from the csv file
    that is used to run the experiment. 
    """
    input: join(config['pilot_dir'], "{sample}.varscan.txt")
    output: join(config['pilot_dir'], "{sample}.varscan.csv")
    params: metadata = config['samples']['file']
    conda: '../envs/r.yml'
    script: "../scripts/add_metadata.R"


rule pilot_aggregate_variants:
    """
    This rule aggregates all of the variants. 
    """
    input: expand(join(config['pilot_dir'], "{sample}.varscan.csv"), sample = ['SSPE_1', 'SSPE_2'])
    output: join(config['pilot_dir'], "pilot_variants.csv")
    run: aggregate_csv(input, output)
