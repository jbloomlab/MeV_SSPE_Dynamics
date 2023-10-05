# Acquiring brain tropism: the spatial dynamics and evolution of a measles virus collective infectious unit that drove lethal subacute sclerosing panencephalitis 

Iris Yousaf, William W. Hannon, Ryan C. Donohue, Christian K. Pfaller, Kalpana Yadav, Ryan J. Dikdan, Sanjay Tyagi, Declan C. Schroeder, Wun-Ju Shieh, Paul A. Rota, Alison F. Feder, Roberto Cattaneo

It is increasingly appreciated that pathogens can spread as infectious units constituted by multiple, genetically diverse genomes, also called collective infectious units or genome collectives. However, genetic characterization of the spatial dynamics of collective infectious units in animal hosts is demanding, and it is rarely feasible in humans. Measles virus (MeV), whose spread in lymphatic tissues and airway epithelia relies on collective infectious units, can, in rare cases, cause subacute sclerosing panencephalitis (SSPE), a lethal human brain disease. In different SSPE cases, MeV acquisition of brain tropism has been attributed to mutations affecting either the fusion or the matrix protein, or both, but the overarching mechanism driving brain adaptation is not understood. Here we analyzed MeV RNA from several spatially distinct brain regions of an individual who succumbed to SSPE. Surprisingly, we identified two major MeV genome subpopulations present at variable frequencies in all 15 brain specimens examined. Both genome types accumulated mutations like those shown to favor receptor-independent cell-cell spread in other SSPE cases. Most infected cells carried both genome types, suggesting the possibility of genetic complementation. We cannot definitively chart the history of spread of this virus in the brain, but several observations suggest that mutant genomes generated in the frontal cortex moved outwards as a collective and diversified. During diversification, mutations affecting the cytoplasmic tails of both viral envelope proteins emerged and fluctuated in frequency across genetic backgrounds, suggesting convergent and potentially frequency-dependent evolution for modulation of fusogenicity. We propose that a collective infectious unit drove MeV pathogenesis in this brain. Re-examination of published data suggests that similar processes may have occurred in other SSPE cases. Our studies provide a primer for analyses of the evolution of collective infectious units of other pathogens that cause lethal disease in humans. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8412085.svg)](https://doi.org/10.5281/zenodo.8412085)

## Snakemake Pipeline

This project was a collaboration between the [Cattaneo Lab](https://www.mayo.edu/research/faculty/cattaneo-roberto-ph-d/bio-00027692), [Bloom Lab](https://research.fredhutch.org/bloom/en.html), and [Feder Lab](https://federlab.github.io/).

This repository contains the Snakemake pipeline that runs the analysis for the study "Acquiring brain tropism: the spatial dynamics and evolution of a measles virus collective infectious unit that drove lethal subacute sclerosing panencephalitis". 

The pipeline can be run by cloning the repository: 

```bash
git clone git@github.com:jbloomlab/MeV_SSPE_Dynamics.git
cd MeV_SSPE_Dynamics
```

Setting up the conda environment specified in the `environment.yml` file: 

```bash
conda env create -f environment.yml
conda activate MeV
```

And finally, running the snakemake pipeline: 

```bash
snakemake --cores 4
```

or, if you're using the `rhino` computing system at Fred Hutch Cancer Cetner, by submitting the following script to sbatch: 

```bash
sbatch run_analysis.bash
```

## Repo Organization

- [`workflow`](workflow/): Contains the code that implements that [Snakemake pipeline](workflow/Snakefile) and the [notebooks](workflow/notebooks/) where the main analysis happens.
- [`config`](config/): Contains some of the data needed to run the pipeline, including the [custom MeVChiTok reference genome](config/ref/MeVChiTok.fa), [annotations](config/annotations.csv), and a [list of sequencing runs](config/samples.csv).
- [`results`](results/): Contains a subset of the results including most importantly the [filtered variants](results/variants/validated_variants.csv) annotated by their haplotype identity.


