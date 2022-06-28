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
| Organism | The name of the organism. *See note* |
| RefGenome | Reference genome accession, if applicable. *See note* |
| BioProject | If applicable. Otherwise any value is acceptable. |
| fq1 | Optional. Path to read 1 for sample |
| fq2 | Optional. Path to read 2 for sample | 

It is important to note that samples are proccessed together based on their `Organism` and `RefGenome` metadata. Thus if you want all of your samples genotyped together, all samples **must all share the same** `Organism` **name and also must share the same** `RefGenome` **value.** 

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

samples: "sample_sheets/<sample_sheet>.csv" # path to sample metadata CSV 
tmp_dir: "tmp/"  # directory path for a temp dir 
split_by_n: True # set to False to split by chromosome/scaffold; set to True to split on runs of Ns within chromosomes/scaffolds.
remote_reads: False #set to True if your reads are stored in a remote Google Cloud bucket
remote_reads_prefix: "" # name of bucket where reads live if above True.
```

To run this out of the box, we reccomend setting split_by_n to True and sentieon to False. 

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

Output can be found in the `results` folder and includes everal key files:

* `results/{SPECIES_NAME}/{ASSEMBLY_NAME}`

The main output of the pipeline is a single VCF with genotype calls for every individual: 

* `results/{SPECIES_NAME}/{ASSEMBLY_NAME}/{SPECIES_NAME}_{ASSEMBLY_NAME}.final.vcf.gz`

By default, this file contains all SNPs and Indels identified and has the basic GATK filters applied. No filtering has been done on the VCF, so it will include all individuals from the sample sheet and all variants identified. The filters are applied as annotations within the VCF file. 

A very simple example for removing filtered sites and only retaining biallelic SNPs is, e.g.:
```
bcftools view -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | ALT="*" | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
```
