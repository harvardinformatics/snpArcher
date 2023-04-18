# Using snpArcher
## Environment Setup
First, we recommend installing Snakemake in a fresh [Mamba](https://github.com/mamba-org/mamba) environment:
```
mamba create -c conda-forge -c bioconda -n snparcher snakemake
mamba activate snparcher
```
Please see the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for detailed installation instructions.

Next, clone the [snpArcher github repo](https://github.com/harvardinformatics/snpArcher) to your machine:
```
git clone https://github.com/harvardinformatics/snpArcher.git
cd snpArcher
```

## Creating a sample sheet
In order to determine what outputs to create, snpArcher requires sample sheet file. This comma separated file contains the required sample metadata about the your samples in order to run the workflow. At a minimum, the snpArcher pipeline requires that each sample have a unique sample name, a reference genome accession or a path to a fasta file, and a SRA accession, or path to two paired end fastq files. 

Below are all of the accepted fields for a sample sheet:
| Field | Description |
| ---- | -------------|
| BioSample | The name of the sample. |
| LibraryName | LibraryID for sample, **must be unique.** |
| Run | The SRR for the sample, if applicable. If not, must be some **unique** value. |
| refGenome | Reference genome accession, if applicable. *See note* |
| refPath | Path to local reference genome, if applicable. *See note* |
| BioProject | If applicable. Otherwise any value is acceptable. |
| fq1 | Optional. Path to read 1 for sample |
| fq2 | Optional. Path to read 2 for sample |
| SampleType | Optional. Triggers postproccesing module. Accepted values are 'include' or 'exclude' |

```{note}
refGenome is always required. refPath is optional, but when specified, a name for the assembly (in refGenome) must also be included. 
```

It is important to note that samples are proccessed together based on their `refGenome` metadata, so **all BioSamples that share a reference genome will ultimately end up in the same final vcf file.** If you are mapping multiple populations / species to a single reference genome, and want separate VCF files for each population / species, you will need to split your final vcf after the pipeline completes, or run multiiple indpendent sample sheets in different results directories. 

If your reads (and, optionally, your local reference genome) are stored in somewhere seperate of the workflow (e.g.: a scratch disk) then you can specify the path to your reads using the `fq1` and `fq2` fields, and the location of your reference genome fasta (*note: must be uncompressed*) in the `refPath` field. 

## Using data from NCBI SRA
If you'd like to reanalyze an existing NCBI SRA BioProject, please follow these instructions to quickly create a sample sheet.

1. Go to the BioProject overview web page on the SRA.
2. In the subheading 'Project Data' there is a table with the columns 'Resource Name' and 'Number of Links'. Click the link in the 'Number of Links' column for the 'SRA Experiments' row. You will be redirected to a search results page.
3. Near the top of the search results page, click the link 'Send results to Run Selector'
4. On the Run Selector page, you can select/deselect samples you'd like to include/exclude in your sample sheet by using the checkboxes.
5. Once you are done selecting samples, click the 'Metadata' button in the 'Download' column in the table near the middle of the page.
6. This will download a a comma separated file called 'SraRunTable.txt'.
7. Open the file in the editor of your choice, and add a column named 'refGenome'. In this column, enter the reference genome accession you want to use for every row in the sheet.
8. Save the sample sheet, it is now ready to use.

### Using local data
If you already have data on your machine, you can use `workflow/write_samples.py` to quickly create a sample sheet. This script takes three arguments:



