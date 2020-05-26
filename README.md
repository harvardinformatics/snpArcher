This is a snakemake pipeline to go from fastq's to VCFs.

TO DO:

clean up snakemake files, partitioning rules into separate .smk files. Follow other guidelines that pros do.
Label rules 01\_, 02\_, etc.

call programs not in your home directory

use profile instead of cluster.json file.

add vcf sanity checks with vcftools.

Also, we should add more preliminary steps that checks and filters fastq files, since errors here may be carried downstream.

We can also add snpEff, but a database must exist, or gff file provided.

make many of these files temporary, but do this later so that entire pipeline doesn't need to be rerun for testing later steps.


## Snakemake pipeline 1: fastq -> BAM

![](workflowScheme_fastq2bam.png)

This is an intentional stopping point to ensure that all BAMs look alright according to the numbers in the "sumstats.txt" file. E.g. make sure dequencing depths are around what you expect etc.

## Snakemake pipeline 2: fastq -> BAM

![](workflowScheme_bam2vcf.png)
