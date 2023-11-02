# Running snpArcher
## Setup
Please refer to our [setup instructions](./setup.md) to prepare the snpArcher environment and requisite files. 
## Test datasets
To test that your environment is properly setup, you can run a quick test with the following command:
```
snakemake -d .test/ecoli --cores 1 --use-conda
```
If this runs without errors, you are ready to go!
## Using the Dry-run option
Snakemake offers the `--dry-run (-n)` CLI option to perform a dry-run of the workflow to show what jobs would be run. We recommend doing this before executing snpArcher to ensure that the sample sheet was setup correctly, and Snakemake has correctly built the workflow DAG.
## Local Execution
Once you have setup the requisite configuration files and sample sheet, executing snpArcher on your local machine is as simple as running the Snakemake command with the number of cores you would like to use. For example, to use 8 cores you would run:
```
snakemake --cores 8 --use-conda
```

### Optional directory setup 
To maintain organization across many different projects, you may consider creating a new directory for each project you run using snpArcher. This way, each of your project directories will contain the configuration files used for that run. Below is an example directory structure:

```
.
├── snpArcher
├── project_1/
│   ├── config/
│   │   ├── config.yaml
│   │   ├── resources.yaml
│   │   └── samples.csv
│   ├── data
│   └── results
└── project_2/
    ├── config/
    │   ├── config.yaml
    │   ├── resources.yaml
    │   └── samples.csv
    └── data
```

When creating a new directory for an analysis, ensure that you copy the `config` directory from the snpArcher directory to your new directory.

Then, to run snpArcher on `project_2` from our example, we would execute the command:
```
snakemake -s ./snpArcher/workflow/Snakefile -d ./project_2 --cores <num cores to use> --use-conda
```

## Cluster Execution
Snakemake [supports most cluster schedulers](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). Here, we provide documentation for SLURM, however please refer to Snakemake's documentation for further details on using other schedulers.
### SLURM
To execute snpArcher on a SLURM cluster, you will need to use the `--profile` Snakemake option. We have included a profile already in, `profiles/slurm`, you will need to edit the `config.yaml` and `cluster-config.yaml` files in this directory. 
#### Profile Setup
1. `config.yaml` defines Snakemake arguments such as number of cores to use and number of jobs to run simulatenously. These values should be adjusted in accordance with your cluster's guidelines.
2. `cluster-config.yaml` defines information about your cluster such as partiton names, number of nodes to use, and runtime. Please edit these for your cluster.
#### Running the workflow
To submit the main snpArcher job to your cluster you will need to call Snakemake in a SLURM job script. We include an example, `run_pipeline.sh`:
```
#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate snakemake
snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm

```



## Cloud Execution
Like cluster execution, Snakemake [supports a number of cloud providers](https://snakemake.readthedocs.io/en/stable/executing/cloud.html). Here we provide documentation for executing using Snakemake's Google Lifesciences integration. Please refer to Snakemake's documentation for details on using other cloud providers.
### Google Lifesciences
Snakemake's integration with the Google Lifesciences (GLS) API allows you to easily run snpArcher on the Google Cloud Platform (GCP). Using this execution mode allows you to take advantage of hundreds or thousands of GCP virtual machine instances. Snakemake manages deploying instances, running jobs, and deleting instances with finished jobs. Furthermore, you can use preemptible instances which are offered at a large cost discount, but can only run for a maximum of 24 hours. 

We include profiles for GLS using preemptible instances so that you can get up and running quickly.
#### Google Credential Setup
In order to use the Google Lifesciences execution option, you must first setup your Google Cloud credentials. Please refer [here](https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html#credentials) for full details.
#### Data setup
To use this execution mode, you must have a Google Storage bucket with your raw data files. This can be achieved by using Google's web interface, or at the command line using [`gsutil`](https://cloud.google.com/storage/docs/gsutil). For example, if we have some data locally like so:
```
.
└── data/
    ├── raw_reads/
    │   ├── samp_1_R1.fq.gz
    │   ├── samp_1_R2.fq.gz
    │   ├── samp_2_R1.fq.gz
    │   ├── samp_2_R2.fq.gz
    │   ├── samp_3_R1.fq.gz
    │   ├── samp_3_R2.fq.gz
    │   ├── samp_4_R1.fq.gz
    │   ├── samp_4_R2.fq.gz
    │   └── ...
    └── ref_genome/
        └── genome.fa
```

We can copy this data to our bucket like so:
```
gsutil cp -r ./data gs://<bucket-name>
```

```{note}
When using cloud execution, do not include the bucket name in any path fields of the sample sheet, such as fq1, fq2, or refPath
```
```{note}
If you are using data hosted on NCBI, you do not need to upload those data to your bucket, snpArcher will handle this for you. However, you still need to create storage bucket to be used for the workflow.
```

Some users may want to store their raw reads in a separate bucket from where the workflow will store files. To do so, you can specify the remote prefix in `config/config.yaml`.

#### Running the workflow
Once your credentials and data are setup, you can run snpArcher using the included profiles. If you are using the default GATK based workflow select `profiles/gls-gatk`, for Sentieon use `profiles/gls-sentieon`. These profiles are set to use preemptible instances and run a maximum of 150 jobs concurrently. These values can be adjusted based on your needs. 

To run the workflow, execute the following command:
```
snakemake --profile <GLS profile> --default-remote-prefix <bucket name>
```

As the workflow runs, Snakemake will print out logging information to the terminal. Please refer [here](https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html#step-5-debugging) for further details.
