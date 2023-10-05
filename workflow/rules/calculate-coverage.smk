"""
### ======= Calculate read coverage ======= ###
This snakemake file uses samtools to calculate the read coverage for each sample.

Since this is most important for variant calling, I am using the realigned bam files as input.

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""


rule samtools_depth:
    """ 
    Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam")
    output: join(config['coverage_dir'], "{accession}", "{accession}.depth")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -m 0 -a -q {params.score} -g DUP {input.bam} > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_depth:
    """ 
    Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{accession}", "{accession}.depth"), accession=samples)
    output: depth = join(config['coverage_dir'], "merged.depth"),
            header = temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "POS\tDP\tAccession"}}1' {output.header} > {output.depth}
        """


rule calculate_insert_size: 
    """ 
    Calculate the insert size using samtools stats
    """
    input: bam = join(config['realign_dir'], "{accession}", "{accession}.sorted.bam")
    output: join(config['coverage_dir'], "{accession}", "{accession}.stats")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools stats {input.bam} | grep ^IS | cut -f 2- > {output}
        """


rule merge_insert_size:
    """ 
    Merge the samtools stats tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{accession}", "{accession}.stats"), accession=samples)
    output: join(config['coverage_dir'], "merged.stats")
    run:
        results = [] 
        for infile in input:
            # Extract the accession name from the file name
            accession = os.path.basename(infile).split(".")[0]
            # Read the file using pandas
            df = pd.read_csv(infile, sep='\t', header=None, names=["insert_size", "total_pairs", "inward_oriented_pairs", "outward_orient_pairs", "other_orientation_pairs"])

            # Calculate the mean insert size
            total_pairs = df["total_pairs"].sum()
            weighted_insert_sizes = (df["insert_size"] * df["total_pairs"]).sum()
            mean_insert_size = weighted_insert_sizes / total_pairs
            # Calculate the weighted variance and then standard deviation
            weighted_variance = ((df["insert_size"] - mean_insert_size) ** 2 * df["total_pairs"]).sum() / total_pairs
            std_dev = np.sqrt(weighted_variance)

            # Append the result
            results.append(pd.DataFrame({"accession": [accession], "mean_insert_size": [mean_insert_size], "std_dev": [std_dev]}))

        # Concatenate the results
        result_df = pd.concat(results)

        # Calculate overall mean of the mean insert sizes and standard deviations
        overall_mean_insert_size = result_df["mean_insert_size"].mean()
        overall_std_dev = result_df["std_dev"].mean()

        # Append these values to the dataframe
        result_df = result_df.append(pd.DataFrame({"accession": ["Overall"], "mean_insert_size": [overall_mean_insert_size], "std_dev": [overall_std_dev]}), ignore_index=True)

        # Save to output file
        result_df.to_csv(output[0], sep='\t', index=False)