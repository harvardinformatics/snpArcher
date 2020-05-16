This is a snakemake pipeline to go from fastq's to VCFs.

Some notes:
the reference genome must be properly indexed for BWA and GATK beforehand, using bwa's "index", picard's "CreateSequenceDictionary", and samtools "faidx".
