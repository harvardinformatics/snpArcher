# Automated short-read mapping and variant calling

This is a suite of snakemake pipelines that takes short-read fastq files, maps them to a reference genome, and calls variants to produce a VCF file along with summaries using an HPC cluster. These pipelines are split into two modular parts: 
    1. reference-based read mapping to produce BAM files (fastq -> BAM) 
    2. using BAM files to call variants in one of several ways (BAM -> VCF)

The workflows here allow users to either start with raw fastq files **or** with BAM files if they have already mapped their short-read data to a reference genome. If you start with raw fastq files, these workflows will not directly give you a VCF file. You must first use the fastq->BAM workflow to create BAM files, and once this has been done, you may then use the BAM->VCF workflows to call variants. This stopping point will force you to inspect the quality of you BAM files (which we facilitate by computing several informative metrics, see below) before proceeding to the variant calling workflows, which involve computationally-expensive tasks and the submission of many jobs.

The first part of this workflow maps short reads to a reference genome using a single short-read aligner: BWA. However, the second part gives you two options for variant calling: GATK4 or freebayes. You may also use both variant-calling programs if you like, but they are currently two separate workflows that need to be run individually by the user. Using both variant callers may be useful if you want to select only high quality variants detected by multiple programs.

A key feature of the variant calling workflows is that we have designed a simple algorithm to split the reference genome into many smaller subsegments that are processed in parallel. These subsegments are flanked by strings of N's in order to avoid edge effects. Such lists of subsegments already exists for some organisms (e.g. Humans), but here we create them ourselves so that these workflows may be used with any non-model organisms.


## Getting started

First clone this repository to copy all the files into a directory of your choosing, and move incto this directory:
```
git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling
cd shortRead_mapping_variantCalling
```

Witin this directory you should see a file named *config.yaml* that stores many variables for the various workflows. These variables include the location of files (e.g. reference genome, fastq files, etc) so that the workflows know where to find them, along with the suffixes of certain files (e.g. "_1.fastq.gz" for raw read files) that allow the programs to identify samples and the correct files to use. Please navigate to the section at the top that contains variables that need to be changed; notes within this file describe these variables. The variables in sections lower in the config.yaml file do not necessarily need to be altered. One of these variables that *may* need to be changed is "minNmer", which is the minimum length of an Nmer (e.g. string of 200 N's) used to break up the genome into smaller intervals to be processed independently (which dramatically speeds up the workflow). The larger the Nmer, the lower the likelihood a pair of reads maps to either side, which may create edge effects when we only consider sub-chromosomal intervals for variant calling. The appropriate Nmer length may also depend on the assembly, as programs differ in how many intervening N's they insert to unite contigs into scaffolds. However, if larger values of "minNmer" are specified, than the algorithm has fewer places to create intervals.

After updating the config.yaml file, you may now run one of the workflows, which gets submitted as a job that itself submits many jobs (a maximum of 1000, but this may be changed). If you are runing the fastq -> BAM workflow, simply type the following on the command line to submit this workflow as a job:
```
sbatch run_fastq2bam.sh
```

The BAM -> VCF workflow currently contains two different options, GATK4 or freebayes. These may both be found in the `run_bam2vcf.sh` file. However, I run just one of these workflows by commenting out the one I do not want to run (with a "#" sign at the beginning of the line). If I want to run GATK4, I comment out the second line of text (below the #SLURM directives) containing the snakemake file `Snakefile_bam2vcf_fb`, which is the freebayes pipeline. Likewise, if I want to run the freebayes workflow I comment out the first line (again, below the #SLURM directives) that contains the snakemake file `Snakefile_bam2vcf_gatk`. After this, I type the following on the command line to submit one of these workflows as a job:
```
sbatch run_bam2vcf.sh
```

Once the workflow is submitted as a job, it will output intermediate and final files in subdirectories depending on the workflow (e.g. *fastq2bam*, *gatk*, or *freebayes*). It may take a while before the workflow does any actual work or submitting of jobs, as conda takes a bit to build the software environment.

The workflows successfully completed if the final summary file (described below) are in the appropriate directory. For the fastq -> BAM workflow, this corresponds to the `bam_sumstats.txt` file, and for the BAM -> VCF workflow this corresponds to the `Combined_hardFiltered.vcf` file along with the files summarizing the VCF: `SNP_per_interval.txt` and `missing_data_per_ind.txt`.

## Description of output files

## Changing the versions of programs

The versions of the various programs may be found in the YAML files in the `envs/` directory. You may update any programs listed under the 'dependencies' heading, replacing the version number with the latest you can find after searching the [Anaconda cloud](https://anaconda.org/).

### Test Data

There are currently two different test datasets that accompany this workflow. The zebrafinch data consists of reads for 3 individuals that map to a genome with 3 scaffolds (each 200kb in length). The Black head duck data consists of reads for 3 individuals that maps to a genome with a single scaffold that gets split (by Nmers) into subintervals.

## Workflow components
### Part 1: fastq -> BAM

![](docs/workflowSchematic_fastq2bam.png)

As mentioned above, we encourage users to inspect the BAM files at this point. Please find many informative metrics in "bam_sumstats.txt". This file contains the following information:
1. Sample name
2. Fraction of reads that passed fastp filter
3. Number of filtered reads
4. Percent of PCR duplicates
5. Percent of paired-end reads that mapped to the genome with mapping quality greater than 20
6. Percent of aligned bases that have base qualities greater than 20
7. Mean sequencing depth per site
8. Number of bases covered *at least* once
9. Logical test for whether your BAM file is valid and ready for variant calling (according to picard's ValidateSamFile tool). If not, check the appropriate _validate.txt file in the 02_bamSumstats dir. "FALSE" values indicate that variant callers may fail downstream, although not necessarily as this validation step is very fussy.

### Part 2: BAM -> VCF

Calling variants can be a time-consuming task for eukaryotic data, but running independent tasks in parallel can speed this process up. There are two approaches to [parallelization](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather): multi-threading and scatter-gather. Multi-threading involves instructions written into the program itself to use multiple cores on a machine (or node on Cannon computing cluster), while the scatter-gather approach involves running multiple copies of the program simultaneously on indepdendent tasks. For instance, we could give GATK 10 cores for multi-threading, or we could divide the genome up into 10 parts and run GATK with 1 core on each of these parts, gathering the results into a single file at the end. This workflow uses the scatter-gather approach, as multithreading is not fully supported at all GATK steps at the moment. 

The scatter-gather approach requires dividing the genome up into segments (e.g. entire chromosomes/scaffolds or even subintervals within scaffolds), and storing these segments in separate "list" files that we later give to GATK to tell it to only work on that particular interval. Ideally, each list file would represent similar fractions of the genome so that all tasks take a similar amount of time. These list files have been created by the Broad for the human genome, but we have written an algorithm into this pipeline that will create these intervals for non-model organisms. Briefly, this algorithm works by finding all regions of the genome with consecutive N's (Nmers), as it is [recommended](https://gatk.broadinstitute.org/hc/en-us/articles/360036823571-ScatterIntervalsByNs-Picard-) that subintervals within scaffolds are flanked by Nmers to avoid problems with variant calling at the edges of intervals. The user can specify in the config file the minimum Nmer length to decide whether to split a scaffold into smaller subintervals. A minimum length of 500bp is probably fine, and you can likely go down to 100bp since this is the length of [padding used](https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals-) when calling variants within exomes.

At its peak of resource consumption, the bam -> VCF workflow submits up to (# samples)X(# list files) jobs, although one may control the maximum number of jobs snakemake submits on the command line.

![](docs/workflowSchematic_bam2VCF.png)

If something in the workflow fails, check the log file and look for the keyword "Error", which should direct you to the specific tasks that failed. It is very possible that some errors may be fixed by simply rerunning the Snakefile, as temporary hardware issues may cause errors, e.g. the computing cluster not responding which can cause Input/Output errors. Occasionally, one step will fail because a previous step  produced truncated output (I have seen fastp produce corrupted files that get fed to bwa, where the error ultimately occurs that snakemake detects). In these cases, these files must be manually removed for snakemake to reproduce them.

Make sure all programs are updated!! GATK actively changing and bugs being fixed all the time!!

To change the resources each task requests, please see the cluster_config.yml file in the subdirectory profiles/slurm/. However, a few rules have their resources specified within the rule specification (in rules subdir) so that the requested memory can be incremented with each attempt. I was not able to get this feature working while specifying memory for that rule within the cluster_config file.

## TO DO:

fir variables continaing directory, ask if they end in "/" otherwise add this!

- ive tried the following to address the problem below, re. resubmitting with many resouces. It seems resources need to be specified in the rule, with the resources keyword, and multiplied by the special 'attempt' variable. However, if any resources are specified within cluster_config.yml under the default, these always override resources specified in the rule and it doesn't work. Moreover, if I instead use a value obtained from a dict, it also doesn't work. Basically the only way I'm able to get things to work now is if I specify the number directly in the rules file.
- it really seems like if there are any job submission parameters defined in cluster_config.yml, either for the specific rule or __default__, it just uses those and ignores any job-specific resource allocations.
resubmit failed jobs; for genomicsdbImport, resubmit with more mem too 
make failed jobs resubmit tasks with more resources! e.g. genomicsDBImport
also resubmit haplotypecaller jobs that fail for inexplicable reasons, maybe with 30% more memory?
resources:
        mem_mb=lambda wildcards, attempt: attempt * 100
        # with --restart-times 3, attempt will take on values 1 thru 3

update this, putting all the options found in profile/slurm/config.yaml

post VCF stuff: number of SNPs, number of filtered SNPs, SFS,rrelatedness (vcftools), PCA (the low depth version), NJ tree, SNPs per bp for each scaffold (or any metric that indicates regions of the genome look bad).



have some recommendations about requesting memory/time, maybe even running it on a few individuals first to see what bam2gvcf takes. recommend serial_requeue for bam2gvcf since it doesn't take long? this is also the step that submits the most jobs by far, so keeping resource requests light is important. Also keep in mind that a few (3) Gb gets subtracted from what you request, since jave needs a few extra to run things in the background. CHANGE SCHEMATICS TO HAVE SUGESTED QUEUE, RESOURCE ALLOCATION, suggest low-pending-time queue for bam2vcf. SOme jobs submitted to serial_requeue may fail for strange reasons (e.g. "ModuleNotFoundError"), but resubmitting them by restarting the snakemake pipeline should do the trick.

make bam2vcf workflow more flexible; e.g. it assumes bams names sample_dedup.bam, but users may have bams with diff names; picard scatterByNs requires indexed genome and dequence dict, but this gets done in fastq2bam. To see how to do this best, you could practice on other datasets, starting at the BAM stage.

have option for low seq depth that uses particular tools? E.g. relatedness also depends on depth

have more input checks? e.g., if you specify wrong suffix, an obscure error comes up in the snakemake rules

how to make pipeline have less variability across runs? make script that takes vary large scaffolds and creates subintervals?

practice having the snakemake job fail in several ways, e.g. timeout, and provide notes on how to restart the workflow, e.g. using the --unlock command and having snakemake delete incomplete files?


run on other datasets and have other people try to use it

report median depth of BAMs as well, use gzipped VCFs 

output log files to separate dir.

can snakemake check if a file is corrupted? if so, remove so that pipeline can be rerun without having to manually remove the file? Saw this for fastp... produced corrupted files but still completed successfully (no exit status greateer than 1?) such that bwa failed.

make streamlined system for eliminating particular samples that just aren't behaving. Currently, you have to remove them from the fastq dir.



add vcf sanity checks with vcftools.

Also, we should add more preliminary steps that checks and filters fastq files, since errors here may be carried downstream.

We can also add snpEff, but a database must exist, or gff file provided.

make many of these files temporary, but do this later so that entire pipeline doesn't need to be rerun for testing later steps.
