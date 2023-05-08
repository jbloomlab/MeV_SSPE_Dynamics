# Analysis Notebooks

This subdirectory contains the `python` and `R` notebooks that perform the bulk of the analysis for this project. Below is a description of what each notebook contains.

I numbered the notebooks based on their order in the analysis. Although order is important for some notebooks (e.g. an output is used in a downstream notebook), in others it is arbitrary.

## Overview

1. `process-variant-calls.Rmd`

   > The goal of this notebook is to determine the quality of the variant calling data. In this notebook, I filter variants based on quality and combine variants across programs. I return a filtered and joined variant dataframe.

2. `determine-main-genotypes.Rmd`

   > The goal of this notebook is to phase mutations by clustering variants that have strongly correlated frequencies across multiple tissue samples. In this notebook, I identify the 'genotypes' `Genome 1` and `Genome 2` and the subclonal genome `Genome 1-1`. These are clusters of mutations that are present in every tissue.

3. `phase-subclonal-mutations.Rmd`

   > The goal of this notebook is to phase more clusters of variants using the same approach we used to define the 'genotypes' `Genome-1` and `Genome-2`. In this notebook I identify a rough set of subclonal haplotypes (or mutation clusters) on the background of G1 and G2.

4. `assign-haplotype-backgrounds.Rmd`

   > The goal of this notebook is to establish a method to genotype SNPs as either belonging to `Genome 1` or to `Genome 2` using reads that 'bridge' between haplotype SNPs and SNPs in `Genome-1` and `Genome-2`. In this notebook I determine the background of subclonal haplotypes.

5. `validate-haplotype-assignments.Rmd`

   > The goal of this notebook is take the haplotypes that we identified and check for problems like homoplasy or other issues. In this notebook, I produce a final set of filtered and validated haplotypes.

6. `prepare-spruce-input.Rmd`

   > The goal of this notebook is to prepare the data for SPRUCE analysis. In this notebook, I calculate the mean haplotype frequency and variance across all samples for each haplotype. I format this data for SPRUCE.

7. `filter-spruce-trees.Rmd`

   > The goal of this notebook is to finalize the phylogenetic relationship of the haplotypes identified by SPRUCE/MACHINA. In this notebook, I filter the SPRUCE trees by removing trees with edges that are not supported by bridging reads.

8. `investigate-driver-mutations.Rmd`

   > The goal of this notebook is to determine how subclonal 'driver' mutations fit on the tree. I also look at the frequeny of the two Fusion C-terminal tail mutations in the brain.

9. `make-phylogenetic-tree.ipynb`

   > The goal of this notebook is to use iqtree to make a phylogenetic tree of the haplotypes.

10. `visualize-phylogenetic-tree.Rmd`

    > The goal of this notebook is to visualize the iqtree with ggplot. In this notebook, I use ggtree to visualize the phylogenetic tree.

11. `cluster-tissues-spatially.Rmd`
    > The goal of this notebook is to take an unbiased approach to show how similar tissue compartments have similar viral populations. In this notebook, I use PCA on the frequency of SNVs in each tissue.
