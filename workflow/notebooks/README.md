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
