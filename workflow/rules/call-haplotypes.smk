### ======= Identify haplotypes with Clique-SNV ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

rule call_haplotypes:
    input: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam")
    output: sam = temp(join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.sam")),
            fasta = join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.fasta"), 
            json = join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.json")
    params:
        prefix = join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted"),
        outdir = join(config['haplotype_dir'], "{accession}"), 
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
            -t {wildcards.t} \
            -tf {wildcards.tf} \
            -el 700 \
            -sp {wildcards.start} \
            -ep {wildcards.stop} \
            -os {wildcards.start} \
            -oe {wildcards.stop} \
            -threads {threads} \
            -outDir {params.outdir} \
            -log 
        """


rule parse_haplotypes: 
    input: fasta = join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.fasta"), 
           json = join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.json"), 
           genome = join(config['index_dir']['bwa'], 'MeVChiTok-SSPE.fa')
    output: join(config['haplotype_dir'], "{accession}", "{accession}-{t}-{tf}-{start}-{stop}.sorted.csv")
    conda: "../envs/pysam.yml"
    script: "../scripts/parse-haplotype-results.py"


def get_intervals(start, end, length, overlap): 
    """
    Get a list of intervals to call haplotypes over. 
    """
    intervals = list()
    
    i = start
    while i < end: 
        if i == start:
            intervals.append((start, i + length))
            i += length 
        elif i + length > end: 
            intervals.append((i - overlap, end))
            i = end
        else:
            intervals.append((i-overlap, i + length-overlap))
            i += length - overlap
    
    return intervals 

# Minimum value for O22 observations
t = [10, 100]
# Minimum frequency of O22 relative to total reads 
tf = [0.05]
# Intervals to call haplotypes over
interval = get_intervals(start = 49, end = 15894, length = 500, overlap = 250)

rule aggregate_haplotypes: 
    input: [f"{config['haplotype_dir']}/{accession}/{accession}-{t}-{tf}-{interval[0]}-{interval[1]}.sorted.csv" for accession, t, tf, interval in list(itertools.product(*[samples, t, tf, interval]))]    
    output: join(config['haplotype_dir'], "haplotypes.csv")
    run: aggregate_csv(input, output)
