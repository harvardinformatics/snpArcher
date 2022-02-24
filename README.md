[![Snakemake](https://img.shields.io/badge/snakemake-≥6.13.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# CCGP WGS Workflow

![CCGP logo](/docs/CCGPhorizontalblue.jpeg)

This is the documentation for the [California Conservation Genomics Project](https://www.ccgproject.org/about-us) whole-genome resequencing workflow. It is managed by the [Corbett-Detig lab](https://corbett-lab.github.io/) at UCSC. 

## Overview

The goal of this workflow is to take whole-genome, short-read, sequencing data and generate variant data by aligning reads and calling variants against a reference genome. It is a restructed workflow from the general pipeline snpArcher, with additional features for running on Google Cloud and using the Sentieon engine.     

We use a Snakemake workflow to manage the various steps in this pipeline. Working familiarity with both the GATK pipeline for variant calling and the Snakemake workflow manager will be useful to run this yourself, but below we outline the required steps to get started

![CCGP Workflow](/docs/ccgp_workflow.png)

## Snakemake environment

Steps to set up environment

## Input formatting

The workflow can be run either on fastq files available locally (i.e. on your server) or can be pulled down from SRA. Example input datasets can be found in `/sample_sheets/`

The other file that needs to be updated can be found under `config/config.yml`:

```
##############################
# Variables you need to change
##############################

samples: "sample_sheets/haliotis_c_workflow.csv"            # name of the sample metadata CSV 
tmp_dir: "tmp/"   # directory path for a temp dir 
split_by_n: True    #set to False to split by chromosome/scaffold; set to True to split on runs of Ns within chromosomes/scaffolds.
sentieon: False  #set to True if you want to use sentieon, False if you want GATK
sentieon_lic: "" #set to path of sentieon license
```

To run this out of the box, we reccomend setting split_by_n to True and sentieon to False. 

The sample sheet with all of the samples to be run in the workflow should be placed in the `samples:` row. 

## Run example data

## Options 

e.g. switch between sentieon on scatter-gather, local vs google cloud

## Output

Output can be found in the `results` folder and includes everal key files:

* `results/{SPECIES_NAME}/{ASSEMBLY_NAME}`

The main output of the pipeline is a single VCF with genotype calls for every individual: 

* `results/{SPECIES_NAME}/{ASSEMBLY_NAME}/{SPECIES_NAME}_{ASSEMBLY_NAME}.final.vcf.gz

By default, this file contains all SNPs and Indels identified and has the basic GATK filters applied. No filtering has been done on the VCF, so it will include all individuals from the sample sheet and all variants identified. The filters are applied as annotations within the VCF file. 

A very simple example for removing filtered sites and only retaining biallelic SNPs is, e.g.:
```
bcftools view -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | ALT="*" | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
```


