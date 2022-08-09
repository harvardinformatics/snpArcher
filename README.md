[![Snakemake](https://img.shields.io/badge/snakemake-≥6.13.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# CCGP WGS Workflow

![CCGP logo](/docs/CCGPhorizontalblue.jpeg)

This is the documentation for the [California Conservation Genomics Project](https://www.ccgproject.org/about-us) whole-genome resequencing workflow. It is managed by the [Corbett-Detig lab](https://corbett-lab.github.io/) at UCSC. 

## Overview

The goal of this workflow is to take whole-genome, short-read, sequencing data and generate variant data by aligning reads and calling variants against a reference genome. It is a restructed workflow from the general pipeline snpArcher, with additional features for running on Google Cloud and using the Sentieon engine.     

We use a Snakemake workflow to manage the various steps in this pipeline. Working familiarity with both the GATK pipeline for variant calling and the Snakemake workflow manager will be useful to run this yourself, but below we outline the required steps to get started

![CCGP Workflow](/docs/ccgp_workflow.png)

## Snakemake environment

[Please see Snakemake's installation instructions.](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation)

## Setup
### Sample metadata sheet
The workflow requires a comma seperated metadata sheet to run. The fields of the sheet are outlined below:
| Field | Description |
| ---- | -------------|
| BioSample | The name of the sample. |
| LibraryName | LibraryID for sample, **must be unique.** |
| Run | The SRR for the sample, if applicable. If not, must be some **unique** value. |
| RefGenome | Reference genome accession, if applicable. Mutually exclusive with refPath. |
| BioProject | Optional. |
| fq1 | Optional. Path to read 1 for sample |
| fq2 | Optional. Path to read 2 for sample |
| refPath | Optional. Path to reference genome. Mutually exclusive with refGenome. |

It is important to note that samples are joint genotyped together based on `refGenome` value.

If your reads are stored in somewhere seperate of the workflow (e.g.: a scratch disk) then you can specify the path to your reads using the `fq1` and `fq2` fields. 

A python script `workflow/write_samples.py` is included to help write the sample sheet for you. This script has four required arguments:
|Argument| Description|
| ------ | ---------- | 
| `-s / --sample_list` | Path to a sample list. One sample per line |
| `-f / --fastq_dir` | Path to directory containing ALL fastq files. It is assumed that each fastq file will contain the sample name uniquely. |
| `-r / --ref` | Path to reference fasta. |
| `-o / --org` | The organism name. |




### Workflow configuration
The other file that needs to be updated can be found under `config/config.yml`:

```
##############################
# Variables you need to change
##############################

samples: "config/ecoli_samples.csv"  # path to the sample metadata CSV 
intervals: True    #Set to True if you want to perform variant calling using interval (split by ns) approach. 
sentieon: False  #set to True if you want to use sentieon, False if you want GATK
sentieon_lic: ".lic" #set to path of sentieon license
remote_reads: False # set if you want reads to be on google cloud storage remote
remote_reads_prefix: "" # set to google bucket name where reads live
```

To run this out of the box, we reccomend setting `intervals` to True and `sentieon` to False. 

The sample sheet with all of the samples to be run in the workflow should be placed in the `samples:` row.

### Run the workflow
Execute the Snakemake workflow by running the command:
```snakemake --use-conda --cores <# of cores to use>```

## Run example data
To run the test data, activate your snakemake conda environment and execute the following command:
`snakemake -d .test/ecoli --use-conda --cores <# of cores to use>`
## Options 

e.g. switch between sentieon on scatter-gather, local vs google cloud, should include bits about slurm here probably.

## Output

Output can be found in the `results` folder and includes several key files:

* `results/{refGenome}/`

The main output of the pipeline is a single VCF with genotype calls for every individual: 

* `results/{refGenome}/final.vcf.gz`

By default, this file contains all SNPs and Indels identified and has the basic GATK filters applied. No filtering has been done on the VCF, so it will include all individuals from the sample sheet and all variants identified. The filters are applied as annotations within the VCF file. 

A very simple example for removing filtered sites and only retaining biallelic SNPs is, e.g.:
```
bcftools view -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | ALT="*" | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
```


