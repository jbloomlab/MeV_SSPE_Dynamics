### ======= Align the reads ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

def get_readgroup(wildcards):
    """ 
    Function to get a string to use as the readgroup for an alignment. 
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    readgroup = samples_df.loc[samples_df.Run == wildcards.accession, ["Sample"]].values.flatten().tolist()[0]
    # Return the readgroup
    return str(readgroup)


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


rule bwa_realign:
    """ 
    Perform short read alignment of the filtered 
    viral readas with `bwa-mem`. Realigning to SSPE consensus.
    
    Sort the aligned reads with samtools sort.
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
    
