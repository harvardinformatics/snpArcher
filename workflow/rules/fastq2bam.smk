rule bwa_map:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        r1 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_1.fastq.gz",
        r2 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_2.fastq.gz",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
    output: 
        bam = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam"),
        bai = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam.bai"),
    params:
        rg = get_read_group
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['bwa_map']['mem']
    log:
        "logs/{refGenome}/bwa_mem/{sample}/{run}.txt"
    benchmark:
        "benchmarks/{refGenome}/bwa_mem/{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam} - && samtools index {output.bam} {output.bai}"

rule merge_bams:
    input:
        merge_bams_input
    output:
        bam = temp("results/{refGenome}/bams/postMerge/{sample}.bam"),
        bai = temp("results/{refGenome}/bams/postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_bams/{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['merge_bams']['mem']
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam} > {log}"

rule dedup:
    input:
        unpack(dedup_input)
    output:
        dedupBam = "results/{refGenome}/bams/{sample}_final.bam",
        dedupBai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    conda:
        "../envs/sambamba.yml"
    resources:
        threads = resources['dedup']['threads'],
        mem_mb = lambda wildcards, attempt: attempt * resources['dedup']['mem']
    log:
        "logs/{refGenome}/sambamba_dedup/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/sambamba_dedup/{sample}.txt"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"