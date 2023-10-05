"""
### ======= Process variants and identify haplotypes ======= ###
This snakemake file contains the main analysis of the project. 
Here, I filter variants, identify haplotypes, and make the paper figures.

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""

rule process_variant_calls:
    """
    The goal of this notebook is to determine the quality of the variant calling data.

    In this notebook, I filter variants based on quality and combine variants across programs.
    """
    input:
        join(config['variant_dir'], "variants.csv"),
        join(config['coverage_dir'], "merged.depth")
    output:
        join(config['notebook_dir'], "process-variant-calls.html")
    params: 
        outcsv=join(config['variant_dir'], "filtered_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/process-variant-calls.Rmd"


rule determine_main_genotypes:
    """
    The goal of this notebook is to phase mutations by clustering variants
    that have strongly correlated frequencies across multiple tissue samples. 

    In this notebook, I identify the 'genotypes' `Genome 1` and `Genome 2`.
    """
    input:
        join(config['pilot_dir'], "pilot_variants.csv"),
        join(config['notebook_dir'], "process-variant-calls.html")
    output:
        join(config['notebook_dir'], "determine-main-genotypes.html")
    params: 
        incsv=join(config['variant_dir'], "filtered_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
        outcsv=join(config['variant_dir'], "genotyped_variants.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/determine-main-genotypes.Rmd"


rule phase_subclonal_mutations:
    """
    The goal of this notebook is to phase more clusters of variants 
    using the same approach we used to define the 'genotypes' `Genome-1` and `Genome-2`.

    In this notebook I identify a rough set of subclonal haplotypes on G1 and G2.
    """
    input:
        join(config['notebook_dir'], "determine-main-genotypes.html")
    output:
        join(config['notebook_dir'], "phase-subclonal-mutations.html")
    params: 
        incsv=join(config['variant_dir'], "genotyped_variants.csv"),
        ancestral=join(config['ref_dir'], "annotated_SSPE_consensus_snps.csv"),
        annotations=config['MeVChiTok']['annotations'],
        outcsv=join(config['variant_dir'], "clustered_variants.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/phase-subclonal-mutations.Rmd"


rule assign_haplotype_backgrounds:
    """
    The goal of this notebook is to establish a method to genotype SNPs
    as either belonging to `Genome 1` or to `Genome 2` using reads that
    'bridge' between haplotype SNPs and SNPs in `Genome-1` and `Genome-2`   

    In this notebook I determine the background of subclonal haplotypes.
    """
    input:
        join(config['bridging_dir'], "bridging_reads.csv"),
        join(config['notebook_dir'], "phase-subclonal-mutations.html")
    output:
        join(config['notebook_dir'], "assign-haplotype-backgrounds.html")
    params:
        incsv=join(config['variant_dir'], "clustered_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
        outcsv=join(config['variant_dir'], "assigned_variants.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/assign-haplotype-backgrounds.Rmd"


rule validate_haplotype_assignments:
    """
    The goal of this notebook is take the haplotypes that we identified
    and check for problems like homoplasy or other issues.

    In this notebook, I produce a final set of filtered and validated haplotypes.
    """
    input:
        join(config['bridging_dir'], "bridging_reads.csv"),
        join(config['notebook_dir'], "assign-haplotype-backgrounds.html")
    output:
        join(config['notebook_dir'], "validate-haplotype-assignments.html")
    params: 
        incsv=join(config['variant_dir'], "assigned_variants.csv"),
        rawvariants=join(config['variant_dir'], "variants.csv"),
        ancestral=join(config['ref_dir'], "annotated_SSPE_consensus_snps.csv"),
        annotations=config['MeVChiTok']['annotations'],
        outcsv=join(config['variant_dir'], "validated_variants.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/validate-haplotype-assignments.Rmd"


rule prepare_spruce_input:
    """
    The goal of this notebook is to prepare the data for SPRUCE analysis. 

    In this notebook, I calculate the mean haplotype frequency and variance 
    across all samples for each haplotype. I format this data for SPRUCE.
    """
    input:
        join(config['coverage_dir'], "merged.depth"),
        join(config['notebook_dir'], "validate-haplotype-assignments.html")
    output:
        join(config['notebook_dir'], "prepare-spruce-input.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        rawvariants=join(config['variant_dir'], "variants.csv"),
        spruce=join(config['variant_dir'], "spruce_input.tsv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/prepare-spruce-input.Rmd"


rule make_SPRUCE_trees: 
    """
    Run SPRUCE on the haplotype frequencies determine a set
    of plausible trees based on the haplotype frequencies.
    """
    input: 
        join(config['notebook_dir'], "prepare-spruce-input.html")
    output:
        join(config['variant_dir'], "spruce_output.tsv")
    params:
        spruce_input=join(config['variant_dir'], "spruce_input.tsv")
    conda: "../envs/spruce.yml"
    shell: 
        """
        generatemutationtrees {params.spruce_input} > {output}
        """


rule filter_spruce_trees:
    """
    The goal of this notebook is to finalize the phylogenetic relationship 
    of the haplotypes identified by SPRUCE/MACHINA.

    In this notebook, I filter the SPRUCE trees by removing trees with 
    edges that are not supported by bridging reads. 
    """
    input: 
        join(config['variant_dir'], "spruce_output.tsv"),
        join(config['bridging_dir'], "bridging_reads.csv")
    output:
        join(config['notebook_dir'], "filter-spruce-trees.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        final_tree=join(config['phylogeny_dir'], "filtered_spruce_tree.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/filter-spruce-trees.Rmd"


rule investigate_driver_mutations:
    """
    The goal of this notebook is to determine how
    subclonal 'driver' mutations fit on the tree. 
    """
    input: 
        join(config['bridging_dir'], "bridging_reads.csv"),
        join(config['notebook_dir'], "filter-spruce-trees.html")
    output:
        join(config['notebook_dir'], "investigate-driver-mutations.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        rawvariants=join(config['variant_dir'], "variants.csv"),
        drivers=config['MeVChiTok']['drivers'],
        annotations=config['MeVChiTok']['annotations'],
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/investigate-driver-mutations.Rmd"


rule make_phylogenetic_tree:
    """
    The goal of this notebook is to use iqtree to make 
    a phylogenetic tree of the haplotypes.
    """
    input: 
        notebook = join(config['notebook_dir'], "investigate-driver-mutations.html"),
        reference = join(config['ref_dir'], "MeVChiTok-SSPE.fa"),
    output:
        alignment = join(config['phylogeny_dir'], "haplotype-sequences.fa"),
        tree = join(config['phylogeny_dir'], "haplotype-sequences.fa.treefile")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        spruce=join(config['phylogeny_dir'], "filtered_spruce_tree.csv")
    conda: "../envs/iqtree.yml"
    notebook: "../notebooks/make-phylogenetic-tree.ipynb"


rule visualize_phylogenetic_tree:
    """
    The goal of this notebook is to visualize the iqtree with ggplot.

    In this notebook, I use ggtree to visualize the phylogenetic tree.
    """
    input:
        tree=join(config['phylogeny_dir'], "haplotype-sequences.fa.treefile")
    output: 
        join(config['notebook_dir'], "visualize-phylogenetic-tree.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        figures=join(config['figure_dir'])
    conda: "../envs/r.yml"
    script: "../notebooks/visualize-phylogenetic-tree.Rmd"


rule cluster_tissues_spatially:
    """
    The goal of this notebook is to take an unbiased approach
    to show how similar tissue compartments have similar viral populations.

    In this notebook, I use PCA on the frequency of SNVs in each tissue.
    """
    input: 
        join(config['notebook_dir'], "visualize-phylogenetic-tree.html")
    output: 
        join(config['notebook_dir'], "cluster-tissues-spatially.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        brain_coords=config['MeVChiTok']['brain_coords'],
        figures=join(config['figure_dir'])    
    conda: "../envs/r.yml"
    script: "../notebooks/cluster-tissues-spatially.Rmd"


rule analyze_strand_origin:
    """
    The goal of this notebook is to see if the strand 
    sense (+ or -) impacts variant calling and coverage.
    """
    input: 
        join(config['notebook_dir'], "visualize-phylogenetic-tree.html"),
        join(config['split_dir'], "merged.split.depth"),
        join(config['split_dir'], "split_read_variants.csv")
    output: 
        join(config['notebook_dir'], "analyze-strand-origin.html")
    params: 
        incsv=join(config['variant_dir'], "validated_variants.csv"),
        annotations=config['MeVChiTok']['annotations'],
        figures=join(config['figure_dir'])    
    conda: "../envs/r.yml"
    script: "../notebooks/analyze-strand-origin.Rmd"


rule clique_snv_analysis:
    """
    The goal of this notebook is to see if G1/G1, G2/G2, and G1/G2 SNV pairs
    are linked for forbidden.
    """
    input: 
        join(config['notebook_dir'], "analyze-strand-origin.html"),
        join(config['bridging_dir'], "bridging_reads.csv")
    output: 
        join(config['notebook_dir'], "clique-snv-analysis.html")
    params: 
        incsv=join(config['variant_dir'], "clustered_variants.csv"),
        figures=join(config['figure_dir'])    
    conda: "../envs/r.yml"
    script: "../notebooks/clique-snv-analysis.Rmd"
