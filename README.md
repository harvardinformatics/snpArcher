## Overview

This fork of the pipeline snpArcher only differs in the following ways:

* custom config decleration in make snakefile for running multiple projects on google cloud simultaneously 
* comment out mappability and coverage filtering code as this does not run on the cloud
* custom snakemake rule for uploading pipeline results to hosting sites to share with CCGP project users
* custom python script for generating CCGP-specific snakemake command that is based on project ID

## Snakemake environment

[Please see Snakemake's installation instructions.](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation)

We recommend installing snakemake into a conda environment called `snakemake`.

## Setup

### Sample metadata sheet

The workflow requires a comma seperated metadata sheet to run. The fields of the sheet are outlined below:
| Field | Description |
| ---- | -------------|
| BioSample | The name of the sample. |
| LibraryName | LibraryID for sample, **must be unique.** |
| Run | The SRR for the sample, if applicable. If not, must be some **unique** value. |
| Organism | The name of the organism. |
| refGenome | Reference genome accession, if applicable. *See note* |
| refPath | Path to local reference genome, if applicable. *See note* |
| BioProject | If applicable. Otherwise any value is acceptable. |
| fq1 | Optional. Path to read 1 for sample |
| fq2 | Optional. Path to read 2 for sample |

*Note:* refGenome is always required. refPath is optional, but when specified, a name for the assembly (in refGenome) must also be included. 

It is important to note that samples are proccessed together based on their `refGenome` metadata, so **all BioSamples that share a reference genome will ultimately end up in the same final vcf file.** If you are mapping multiple populations / species to a single reference genome, and want separate VCF files for each population / species, you will need to split your final vcf after the pipeline completes, or run multiiple indpendent sample sheets in different results directories. 

If your reads (and, optionally, your local reference genome) are stored in somewhere seperate of the workflow (e.g.: a scratch disk) then you can specify the path to your reads using the `fq1` and `fq2` fields, and the location of your reference genome fasta (*note: must be uncompressed*) in the `refPath` field. 

A python script `workflow/write_samples.py` is included to help write the sample sheet for you. This script has three required arguments:
|Argument| Description|
| ------ | ---------- |
| `-s / --sample_list` | Path to a sample list. One sample per line |
| `-f / --fastq_dir` | Path to directory containing ALL fastq files. It is assumed that each fastq file will contain the sample name uniquely. |
| `-o / --org` | The organism name. |

Additionally, one of the following options must also be specified:
|Argument| Description|
| ------ | ---------- |
| `-r / --ref` | Path to reference fasta. |
| `-a / --acc` | NCBI accession of reference. |

### Workflow configuration

The other file that needs to be updated can be found under `config/config.yml`:

```
##############################
# Variables you need to change
##############################

samples: "config/samples.csv" # path to the sample metadata CSV
resource_config: "config/resources.yaml" # path to resources yaml config
final_prefix: "" # prefix for final output files
intervals: True #Set to True if you want to perform variant calling using interval approach.
sentieon: False #set to True if you want to use sentieon, False if you want GATK
sentieon_lic: "" #set to path of sentieon license
remote_reads: False # Set True if reads are in a Google Bucket seperate from --default-remote-prefix.
remote_reads_prefix: "" # set to google bucket prefix where reads live
bigtmp: "" #Set to a path with lots of free space to use for commands that require large amounts of temp space; defaults to system tmpdir if empty
cov_filter: True #set to True if you want to include coverage thresholds in the callable sites bed file (default uses mappability only)

```
You are required to set a `final_prefix`: your final vcf file will be called {final_prefix}_final.vcf. To run out of the box, we recommend leaving `intervals` True and `sention` and `remote_reads` (for cloud operations) False. You will need to change `samples` to point to your sample sheet. You may need to set `bigtmp` to something depending on your system configuration, as some steps can requires 100+ Gb of temp space for large datasets; `tmp/` is a decent option, which will create a temp directory in the directory you launch snpArcher from. 

The default `cov_filter` is very loose, removing only regions of the genome with 0 coverage and excessively high coverage (>10000x). You can modify these parameters in the config, as needed. The default options are also optimized for low coverage (<10x) data. If you have >15x coverage, you probably want to use the high coverage parameters instead, which can also be changed in the config. 

Increasing the `num_gvcf_intervals` will scale the pipeline wider, so there will be more, shorter jobs. Decreasing `num_gvcf_intervals` will create fewer, longer jobs. The optimal setup will depend on your HPC system. 

### Run the workflow

Execute the Snakemake workflow by running the command:
`snakemake --use-conda --cores <# of cores to use>`

If you are running on an HPC system, you'll need to include a profile for, e.g. `slurm` or `sge`:
`snakemake --profile profiles/slurm`

We include a slurm profile with snpArcher, but you will need to modify the queue names, at a minimum, for your institution. We welcome contributions of other profiles. The `run_pipeline.sh` script is an example of how you might launch a snakemake job using slurm, and assuming you have snakemake installed in a conda environment called snakemake, although you will need to modify the partition names and possibly the resource requests depending on your local environment (and we would recommend changing the out and err files to something more informative, e.g. {final_prefix}-sm-%j.out and {final_prefix}-sm-%j.err).

## Run example data

To run the test data, activate your snakemake conda environment and execute the following command:
`snakemake -d .test/ecoli --use-conda --cores <# of cores to use>`

## Output

Output can be found in the `results` folder and includes everal key files:

- `results/{ASSEMBLY_NAME}`

The main output of the pipeline is a single VCF with genotype calls for every individual:

- `results/{ASSEMBLY_NAME}/{final_prefix}_final.vcf.gz`

By default, this file contains all SNPs and Indels identified and has the basic GATK filters applied. No filtering has been done on the VCF, so it will include all individuals from the sample sheet and all variants identified. The filters are applied as annotations within the VCF file.

A very simple example for removing filtered sites and only retaining biallelic SNPs is, e.g.:

```
bcftools view -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | ALT="*" | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
```
