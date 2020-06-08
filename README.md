This is a snakemake pipeline to go from fastq's to VCFs in 2 modular parts.

The first part uses (unfiltered) raw fastq files to produce BAM files, which may be used with a variety of downstream variant callers. Users are encouraged to inspect these BAM files before proceding so that time (human or machine) is not wasted calling variants on defective BAM files (e.g. sample contaminatation, reference sequence too diverged, raw reads of unexpectedly low quality, PCR step from the DNA library preparation caused extreme biases, etc.).

Currently, the second part of our workflow uses GATK4 to call variants, although one may take their BAMs from part one and use other programs. In the future, we may create different versions of this second part to support a variety of popular variant callers so that users may, for example, take the intersection of the results from multiple programs.



The following is a schematic diagram of the workflow:

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

Calling variants can be a time-consuming task for eukaryotic data, but running independent tasks in parallel can speed this process up. There are two approaches to [parallelization](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather): multi-threading and scatter-gather. Multi-threading involves instructions written into the program itself to use multiple cores on a machine (or node on Cannon computing cluster), while the scatter-gather approach involves running multiple copies of the program simultaneously on indepdendent tasks. For instance, we could give GATK 10 cores for multi-threading, or we could divide the genome up into 10 parts and run GATK with 1 core on each of these parts, gathering the results into a single file at the end. This workflow uses both approaches, since the multi-thread option of GATK by itself (currently only supported by the HaplotypeCaller step) has dminishing returns (*verify this*). Thus, we also use the scatter-gather approach, which allows us to also speed up the process by parallelizing parts of GATK that do not support multithreading.

The scatter-gather approach requires dividing the genome up into segments (e.g. chromosomes or scaffolds), and storing these segments in separate "list files" that we later give to GATK. Ideally, each list file would represent similar fractions of the genome. For instance, one list file may contain only the name of one large scaffold, while another may contain two shorter scaffolds that have a cumulative length similar to the larger scaffold in the previous file. To facilitate this approach, we provide a script (createListsFromFasta.py) that takes as input a fasta file and outputs many list files. The user gives two parameters to this script: the number of target segments to break the genome into and the maximum number of segments per list file. This second parameter is specified because incomplete assemblies may have thousands of smaller scaffolds, such that breaking the genome into, say, 20 segments may put thousands of small scaffolds into a particular list file. We have some experience suggesting that thousands of scaffolds in a list file may slow down particular steps of the GATK workflow.


At its peak of resource consumption, the bam -> VCF workflow submits up to (# samples)X(# list files) jobs, although one may control the maximum number of jobs snakemake submits on the command line.

![](docs/workflowSchematic_bam2VCF.png)

If something in the workflow fails, check the log file and look for the keyword "Error", which should direct you to the specific tasks that failed. It is very possible that some errors may be fixed by simply rerunning the Snakefile, as temporary hardware issues may cause errors, e.g. the computing cluster not responding which can cause Input/Output errors. Occasionally, one step will fail because a previous step  produced truncated output (I have seen fastp produce corrupted files that get fed to bwa, where the error ultimately occurs that snakemake detects). In these cases, these files must be manually removed for snakemake to reproduce them.

## TO DO:

output log files to separate dir.

can snakemake check if a file is corrupted? if so, remove so that pipeline can be rerun without having to manually remove the file? Saw this for fastp... produced corrupted files but still completed successfully (no exit status greateer than 1?) such that bwa failed.

make streamlined system for eliminating particular samples that just aren't behaving. Currently, you have to remove them from the fastq dir.

post VCF stuff: relatedness (vcftools), missingness (vcftools), PCA, NJ tree

use profile instead of cluster.json file.

add vcf sanity checks with vcftools.

Also, we should add more preliminary steps that checks and filters fastq files, since errors here may be carried downstream.

We can also add snpEff, but a database must exist, or gff file provided.

make many of these files temporary, but do this later so that entire pipeline doesn't need to be rerun for testing later steps.
