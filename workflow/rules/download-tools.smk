### ======= Download and assembly of tools (independent of `conda`) ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule get_varscan:
    """ 
    Download Varscan into Tools directory. wget from github. 
    """
    output: join(config['tools'], "VarScan.v2.4.0.jar")
    params: http=config['varscan']
    shell: "wget -O {output} {params.http}"


rule get_SnpEff:
    """ 
    Download and build SnpEff with annotations from gtf file. 
    """
    output: join(config['tools'], 'snpEff/snpEff.jar')
    params:
        installdir=config['tools'],
        installpath=config['snpeff'],
        tmp=join(config['tools'], 'snpEff_latest_core.zip'),
        datadir=join(config['tools'], 'snpEff/data')
    shell:
        """
        # Download SnpEff
        wget {params.installpath} -O {params.tmp}

        # Unzip and install
        unzip {params.tmp} \
            -d {params.installdir} \
            && rm -rf {params.tmp}

        # Set path variable to snpEff dir
        snpeff={output}

        # Make the data directory to install viral genomes
        mkdir -p {params.datadir}
        """



rule build_SnpEff:
    """ Build the annotation repository for the genome used to call variants. 
    """
    input: 
        snpeff=join(config['tools'], 'snpEff/snpEff.jar'),
        ref=join(config['ref_dir'], 'MeVChiTok-SSPE.fa'),
        gff=join(config['gff_dir'], 'MeVChiTok.gff')
    output: 
        virusdir=directory(join(config['tools'], 'snpEff/data/MeVChiTok'))
    params:
        contig="MeVChiTok",
        config=join(config['tools'], 'snpEff/snpEff.config')
    conda: '../envs/java.yml'
    shell:
        """
        # Make the virus specific dir
        mkdir -p {output.virusdir}

        # Move the genome fasta file
        cp {input.ref} {output.virusdir}/sequences.fa

        # Move the gtf file
        cp {input.gff} {output.virusdir}/genes.gff

        # Add the genome build to the end of the file
        echo "{params.contig}.genome: {params.contig}" >> {params.config}

        # Build the new database
        java -Xss100M -Xmx8g -jar {input.snpeff} build -gff3 -v {params.contig}
        """

