# Analysis Notebooks

This subdirectory contains the `python` and `R` notebooks that analyze the sequencing data and variant calls from the pipeline. Below is a description of the analysis performed in each notebook.

The notebooks are numbered based on their order in the analysis. While the order is important for some notebooks (e.g. an output is used in a downstream notebook), in others it is arbitrary.

## Overview

1. `quality-control.Rmd`

> In this notebook, I assess the quailty of variant calls by comparing and combining the variants called by different tools (`lofreq` and `varscan`). I also make **extended data figure 3** of the paper where I plot the location of mutations in each genome highlighting regions with low coverage.

2. `identify-backgrounds.Rmd`

> In this notebook, I make a first pass at identifying variants that are on the same viral haplotype by clustering mutations by their correlation in frequency across _every tissue sample_. Here, I identify `Genome 1` and `Genome 2`. I also identify the `Genome 1` sub-haplotype that's present in every sample but the Frontal Cortex 2 specimen.

3. `cluster-subclonal-mutations.Rmd`

> In this notebook, I take the analysis in the previous notebook farther to see if we can cluster mutations by frequency even if they aren't found in every tissue sample. I use this analysis to identify subclonal haplotypes that are on the background of either `Genome-1` or `Genome-2`. I also assign mutations that are missing from the original `Genome-1` or `Genome-2` assignment due to low depth in one or more tissue. Finally, I made `figure 2`, `figure 4a`, and `figure 4b`.

4. `genotype-subclonal-snps.Rmd`

> In this notebook, I use bridging reads and a probabalistic approach written by Alison Feder to assign SNPs to either `Genome 1` or `Genome 2`. This can be used to determine the most likely background for the sub-haplotypes identified in the previous notebook. _TODO: we need to refine the bridging reads approach_

5. `validate-subclonal-haplotypes.Rmd`

> In this notebook, I validate and organize the nomenclature for the subclonal haplotypes that were identified in the previous notebook `cluster-subclonal-mutations.Rmd` and assigned to `Genome 1` or `Genome 2` in the previous notebook `genotype-subclonal-snps.Rmd`. I used bridging reads over any given pair of SNPs to do this validation along with checking for violations in our assumptions about the sum of allele frequencies in each tissue specimen. I also make `extended figure 4` and `figure 5a`.

6. `prepare-spruce-input.Rmd`

> In this notebook, I prepare the reads for `SPRUCE` as implemented by the software program `MACHINA`. I try different approaches to add 'uncertainty' around our estimates for the mean haplotype frequecny in each tissue. I also generate `figure 5b` and `figure 5c`.

7. `analyze-spruce-output.Rmd`

> In this notebook, I parse the output of `SPRUCE` to examine the trees that fit our haplotype data. Then, I attempt to use bridging reads to filter down this list of trees into those that are possible given our data.

---

6. `spatial-tissue-clustering.Rmd`

> In this notebook, I examined how variant populations in different tissues relate to one another in the context of space. For example, are nearby tissues more similar than far away tissues?

7. `mutation-effect-distribution.Rmd`

> In this notebook, I look for patterns in the distribution of mutations in different haplotypes. For example, are there clear clusters of ADAR mutations?

8. `haplotype-phylogenetic-tree.Rmd`

> In this notebook, I work on different ways to visualize the phylogenetic tree of haplotypes.

9. `tissue-level-diversity.Rmd`

> In this notebook, I calculate nucleotide diversity in different tissue compartments.

10. `polish-current-haplotypes.Rmd`

> In this notebook, I use bridging reads to see if I can assign unclustered mutations to our current haplotypes.
