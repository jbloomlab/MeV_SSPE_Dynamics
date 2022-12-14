# Measles SSPE Spatial Dynamics

This repository contains an analysis of within host viral haplotypes of Measles with the aim of constructing a spatial
map of measles diversification in the brain.

## Methods

We performed all bioinformatic analysis of the RNA sequencing of tissue samples isolated from the brain using a computational pipeline created with the [Snakemake workflow language](https://snakemake.readthedocs.io/en/stable/). This pipeline can be found at https://github.com/jbloomlab/MeV_SSPE_Dynamics.

Before aligning reads from the 15 brain isolates, we processed reads by trimming adaptor sequences and low quailty bases with [`fastp`](https://github.com/OpenGene/fastp). After trimming the reads, we extracted the viral reads from human reads using a kmer matching approach from [`BBduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

Following preprocessing, we aligned the viral reads to a reference genome created with a composite sequence that includes all the nucleotides represented in 12 isolates of the MeV genotype D – the only genotype known to have circulated in Central America when the patient may have been infected – using [`BWA`](https://bio-bwa.sourceforge.net/). We took this first alignment and identified the nucleotide mutations found in greater than 12/15 of the tissue samples at more than 90% allele frequency. We used these probably 'ancestral' mutations to build a consensus reference genome of the SSPE samples. We then realigned the preprocessed measles reads to this SSPE consensus genome.

Using the realigned reads, we called variants using two variant calling programs – [`lofreq`](https://csb5.github.io/lofreq/) and [`varscan2`](https://varscan.sourceforge.net/). We restricted our analysis to variants present at more than 2% frequency in more than 200 reads. The variants identfied by these two programs were highly concordant. We combined the variants from both programs by filling in variants that were missing in one variant caller if that variant was present in more than one tissue.

Using the filtered and collated variant calls, we set out to determine the linkage between variants. To do this, we first took variants that were present in all 15 tissue samples and created a correlation matrix using the frequency of each variant in every tissue. We then used k-medoids clustering on the pearson correlation coefficients from this matrix to group mutations that co-varied in frequency across the tissue samples. Mutations that co-vary in frequency are likely to be linked on the same viral haplotype.

To vallidate these haplotypes, we counted the occurence of reads that bridged positions with multiple variants in our putitive haplotypes. Using these bridging reads, combined with a frequency-based exclusion approach, we vallidated our haplotypes and determined the phylogenetic relationship between them.

_You can see a graphical represenation of this pipeline [here](dag.pdf)_

## Repo Organization

`config`: Contains the file that configure the analysis (i.e., the reference genome)
`workflow`:

- `envs`: Conda files for creating environments
- `rules`: Snakemake files with pipeline rules
- `script`: Scripts used by the pipeline
- `notebooks`: Notebooks with analysis and figures
  `results`:
- `variants`: Variant calls (unprocessed)[/results/variants/variants.csv], (filtered)[/results/variants/filtered_variants.csv], (genotyped into G1 and G2)[/results/variants/genotyped_variants.csv], and (clustered into haplotypes)[/results/variants/clustered_variants.csv]
- `notebooks`: HTML files with the knit notebook files and notes.
