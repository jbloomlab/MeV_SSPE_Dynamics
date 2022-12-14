# Analysis Notebooks

This subdirectory contains the python and R notebooks that analyze the sequencing data and variant calls from the pipeline. Below is a description of the analysis performed in each notebook.

The notebooks are numbered based on their order in the analysis. While the order is important for some notebooks (e.g. an output is used in a downstream notebook), in others it is arbitrary.

## Overview

1. `quality-control.Rmd`

> In this notebook, I assess the quailty of variant calls by comparing and combining the variants called by different tools (`lofreq` and `varscan`).

2. `identify-backgrounds.Rmd`

> In this notebook, I make a first pass at identifying variants that are on the same viral haplotype by clustering mutations by their correlation in frequency. Here, I identify `Genome 1` and `Genome 2`. I also identify the `Genome 1` subhaplotype `Genome 1-1`.

3. `cluster-subclonal-mutations.Rmd`

> In this notebook, I take the analysis in the previous notebook farther to see if we can cluster mutations by frequency even if they aren't found in every tissue sample. I use this analysis to identify subclonal haplotypes.

4. `genotype-subclonal-snps.Rmd`

> In this notebook, I use bridging reads and a probabalistic approach written by Alison Feder to assign SNPs to either `Genome 1` or `Genome 2`. This also assigns haplotypes to their most like genotype.

5. `establish-haplotype-relationship.Rmd`

> In this notebook, use a combination of the 'pigeon-hole principle' and bridging reads to determine how the clusters of mutations (putitive haplotypes) are likely related to one another.

6. `spatial-tissue-clustering.Rmd`

> In this notebook, examine how variant populations in different tissues relate to one another in the context of space. For example, are nearby tissues more similar than far away tissues?

7. `mutation-effect-distribution.Rmd`

> In this notebook, I look for patterns in the distribution of mutations in different haplotypes. For example, are there clear clusters of ADAR mutations?

8. `haplotype-phylogenetic-tree.Rmd`

> In this notebook, I work on different ways to visualize the phylogenetic tree of haplotypes.

9. `tissue-level-diversity.Rmd`

> In this notebook, I calculate nucleotide diversity in different tissue compartments.

10. `polish-current-haplotypes.Rmd`

> In this notebook, I use bridging reads to see if I can assign unclustered mutations to our current haplotypes.
