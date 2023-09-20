"""This analysis checks whether the sense of the viral RNA (+/-) impacts the frequency of variants and the coverage over the genome."""

## ======== Split the aligned BAM file into both forward and reverse reads. ======== ##  
rule split_bam:
    """
    Split the BAM files into forward and reverse reads, and then sort and index each of them.
    """
    input:
        bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam"),
        bai = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam.bai")
    output:
        bam = join(config['split_dir'], "{accession}", "{accession}.{direction}.sorted.bam"),
        bai = join(config['split_dir'], "{accession}", "{accession}.{direction}.sorted.bam.bai"),
        tmp_bam = temp(join(config['split_dir'], "{accession}", "{accession}.{direction}.tmp.sorted.bam"))
    params:
        flag = lambda wildcards: "99" if wildcards.direction == "fwd" else "147"
    threads: config['threads']['max_cpu']
    conda: '../envs/align.yml'
    shell:
        """
        # Split into forward and reverse reads
        samtools view -b -f {params.flag} {input.bam} > {output.tmp_bam}
        
        # Sort and index the forward reads BAM file
        samtools sort -o {output.bam} {output.tmp_bam}
        samtools index {output.bam} {output.bai}
        """
## ======== Split the aligned BAM file into both forward and reverse reads. ======== ##  



## ======== Calculate the depth over each position for both forward and reverse reads. ======== ##  
rule samtools_split_depth:
    """ 
    Calculate the depth over each position filtering by the phred base score. 
    """
    input: 
        bam = join(config['split_dir'], "{accession}", "{accession}.{direction}.sorted.bam"),
    output: 
        coverage = join(config['split_dir'], "{accession}", "{accession}.{direction}.depth"),
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -m 0 -a -q {params.score} -g DUP {input.bam} > {output.coverage}
        sed -i "s/$/\t{wildcards.accession}\t{wildcards.direction}/" {output.coverage} 
        """


rule merge_split_depth:
    """ 
    Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: 
        expand(join(config['split_dir'], "{accession}", "{accession}.{direction}.depth"), accession=samples, direction=['fwd', 'rev']),
    output: 
        depth = join(config['split_dir'], "merged.split.depth"),
        header = temp(join(config['split_dir'], "merged.split.depth.tmp")),
    shell:
        """
        cat {input} > {output.header}
        awk 'BEGIN{{print "POS\tDP\tAccession\tDirection"}}1' {output.header} > {output.depth}
        """
## ======== Calculate the depth over each position for both forward and reverse reads. ======== ##  

## ======== Call variants on the split reads for both forward and reverse reads. ======== ##
rule split_read_samtools_mpileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.

    This is necessary for variant calling with VarScan.
    """
    input: bam = join(config['split_dir'], "{accession}", "{accession}.{direction}.sorted.bam"),
           genome = join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa')
    output: join(config['split_dir'], "{accession}", "{accession}.{direction}.mpileup.txt")
    params: score = config['BQ'],
            maxdepth = config['maxdepth']
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d {params.maxdepth} -E --excl-flags UNMAP,SECONDARY,QCFAIL -q {params.score} -Q {params.score} -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"


rule split_read_varscan_calling: 
    """ 
    SNP calling with Varscan.
    
    Parameters are controlled from the config file. 

    Except the strand filter, which is set to 0 to remove the filter. 
    """
    input: 
        varscan = join(config['tools'], "VarScan.v2.4.0.jar"),
        pileup = join(config['split_dir'], "{accession}", "{accession}.{direction}.mpileup.txt")
    output: 
        variants = join(config['split_dir'], "{accession}", "{accession}.{direction}.varscan.vcf")
    params:
        minimum_coverage = config['min_coverage'],
        minumum_supporting_reads = config['min_reads_supporting'],
        minimum_base_quality = config['BQ'],
        minimum_variant_freq = config['min_allele_frequency'],
        strand_filter = 0
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


rule split_read_lofreq_calling: 
    """ 
    Call variants (SNP and indels) with the reference.
    
    Remove the defaut filter beacuse of the strand bias filter.
    """
    input: bam = join(config['split_dir'], "{accession}", "{accession}.{direction}.sorted.bam"),
           genome = join(config['index_dir']['samtools'], 'MeVChiTok-SSPE.fa')
    output: join(config['split_dir'], "{accession}", "{accession}.{direction}.lofreq.vcf")
    params: 
        maxdepth = config['maxdepth']
    conda: '../envs/lofreq.yml'
    threads: config['threads']['max_cpu']
    shell:
        """
        if [[ $(samtools view {input.bam} | head -n 5) ]]; then
            lofreq call-parallel --pp-threads {threads} \
                -f {input.genome} \
                -N \
                -no-default-filter \
                -d {params.maxdepth} \
                {input.bam} \
                -o {output}
        else
            touch {output}
        fi
        """ 


rule annotate_split_read_vcf:
    """
    Use the program SnpEff to annotate the effect of 
    mutations using the custom virus genome as reference.

    This tool is set up in the download_tools.smk file.
    """
    input:
        virusdir = join(config['tools'], 'snpEff/data/MeVChiTok'),
        vcf = join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.vcf")
    output: join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.vcf")
    conda: "../envs/java.yml"
    params:
        snpEff = join(config['tools'], 'snpEff/snpEff.jar'),
        config = join(config['tools'], "snpEff/snpEff.config"),
        contig = "MeVChiTok"
    shell:
        """
        java -jar {params.snpEff} -c {params.config} -noStats \
            -v {params.contig} {input.vcf} > \
            {output}
        """


rule split_read_vcf_to_table:
    """
    Convert the VCF files to tables for easy data
    analysis in R or Python.
    """
    input: join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.vcf")
    output: join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.txt")
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


rule add_split_read_metadata: 
    """
    This rule adds in the metadata from the csv file
    that is used to run the experiment. 
    """
    input: join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.txt")
    output: join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.csv")
    params: metadata = config['samples']['file']
    conda: '../envs/r.yml'
    script: "../scripts/add_metadata.R"


rule aggregate_split_read_variants:
    """
    This rule aggregates all of the variants. 
    """
    input: expand([join(config['split_dir'], "{accession}", "{accession}.{direction}.{caller}.ann.csv")], accession=samples, caller=['lofreq', 'varscan'], direction = ['fwd', 'rev'])
    output: join(config['split_dir'], "split_read_variants.csv")
    run: aggregate_csv(input, output)


