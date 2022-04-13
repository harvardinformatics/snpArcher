# Automated short-read mapping and variant calling

A Snakemake pipeline that wraps GATK and other tools to facilitate variant calling in any organism.

## How to run

### 0.) Install snakemake, if you haven't already
Please follow the [installation via conda instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install snakemake in a conda environment called "snakemake". The script that ultimately runs the snakemake command is preceded by a line with `conda activate snakemake`, so the exact environment name is necessary unless you change this script. If you already have snakemake installed, you can remove this line with `conda activate snakemake` from all the \*sh files.

### 1) Download code
First clone this repository and move into the new directory:
```
git clone https://github.com/harvardinformatics/snpArcher
cd snpArcher
```

### 2) Set values of important variables (config.yaml)
Witin this directory a file named `config.yaml` stores many variables. The most important is the sample sheet, which controls what the pipeline will process. This is a CSV file that contains the sample metadata (following the format of the example data sheet ```samples.csv```). Note that you can use local fastq files by specifying a path instead of an SRA accession, and you can use a local reference genome by adding a `RefPath` column to the sample sheet that gives the path to the reference genome. More complete example sample sheets are coming.

If you are calling data on SRA, make sure that there is one line in the sample sheet for each SRA Run, and that the accessions point to Runs not Experiments or Projects. The code to download data from SRA assumes Run accessions.

### 3) Set the resources to request for various steps (resources.yaml)
The `resources.yaml` file may be changed to increase the amount of requested memory (in Megabytes) or the number of threads for the steps that support multi-threading. Not all steps in the workflows are included here, so these use the default amount of resources. **NOTE**: if any job fails, it gets resubmitted with increased memory calculated as (*attempt number*)\*(initial memory). You may also need to change the values in the `profiles/slurm/cluster_config.yml` to specify appropriate partitions for your HPC environment.

### 4) Are you alright with default number of jobs to submit to run simultaneously?
There's a file in the `profiles/slurm` directory called `config.yaml` which contains various options for the workflow (this setup is from using [profiles](https://github.com/Snakemake-Profiles)). The most important is `jobs` at the top. If your workflow needs to submit ~10k jobs overall and many of them can be run in parallel (e.g. making GVCFs from BAM files for each sample), then this `jobs` variable determines how many jobs the workflow will submit at any given time. The default is 1000, meaning if 1k jobs are sitting in the queue (running or pending), it will not submit more. If you are concerned about your fairshare score decreasing dramatically because of this (e.g. you have 80k jobs to submit overall because you have many samples), set `jobs` to something smaller, such as 300. This will of course make the workflow take longer but will leave resources for your colleagues!

### 5) Submit workflow!
After updating the config.yaml file, you may now run one of the workflows, which gets submitted as a job that itself submits many jobs (max of 1000, see step 4 above if you want to change this). Once the workflow is submitted as a job, it may take a while to build the software environment before it does anything.

Note that the `run_pipeline.sh` script provides an example of what you can submit to run the workflow, while the `run_pipeline_update.sh` gives an example for how to rerun a the workflow with additional samples or other config changes.
