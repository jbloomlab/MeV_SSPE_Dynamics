### ======= Import runs as fastq formatted files ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
#

# Read in the samples dataframe.
SAMPLE_DF = pd.read_csv(config['samples']['file'])

rule fetch_fastq:
    """ 
    Move the fastq files into the results directory. 
    """
    output: join(config['fastq_dir'], "{accession}", "{accession}.fastq.gz")
    params: path = lambda wildcards: SAMPLE_DF.loc[(SAMPLE_DF.Run == wildcards.accession)].Path.item()
    shell: "cp {params.path} {output}"

