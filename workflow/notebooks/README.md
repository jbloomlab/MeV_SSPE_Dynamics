# Analysis Notebooks

This subdirectory contains the `python` and `R` notebooks that perform the bulk of the analysis for this project. Below is a description of what each notebook contains.

I numbered the notebooks based on their order in the analysis. Although order is important for some notebooks (e.g. an output is used in a downstream notebook), in others it is arbitrary.

## Overview

1. `01-process-variant-calls.Rmd`

> In this notebook, I assess the quailty of variant calls by comparing and combining the variants identified by two different variant calling tools (`lofreq` and `varscan`).

2. `02-determine-main-genotypes.Rmd`

> In this notebook, I identify the main viral genotypes by clustering mutations that were identified in **every single tissue sample** at a detectable frequency (>= 2%). This approach identified 3 viral haplotypes/genotypes that were present in every sample; `Genome 1`, `Genome 2`, and a `Genome 1` subhaplotype that's dominant in all samples but the frontal cortex 2 sample.

3. `03-phase-subclonal-mutations.Rmd`

> In this notebook, I use the approach developed in the previous notebook phase viral mutations into haplotypes even if those mutations weren't identified in every tissue sample.

4. `04-assign-haplotype-backgrounds.Rmd`

> In this notebook, I use bridging reads and a probabalistic approach developed by Alison Feder to determine which of the main genotypes each of the haplotypes identified in the previous notebook was descended from.

5. `05-validate-haplotype-assignments.Rmd`

> In this notebook, I validated and refined the haplotypes that I phased in the previous notebook. I also standardize the nomenclature of the haplotypes going forward. 

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
