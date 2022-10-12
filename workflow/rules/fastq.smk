rule get_fastq_pe:
    output:
        temp("results/data/fastq/{refGenome}/{sample}/{run}_1.fastq.gz"),
        temp("results/data/fastq/{refGenome}/{sample}/{run}_2.fastq.gz")
    params:
        outdir = os.path.join(workflow.default_remote_prefix, "results/data/fastq/{refGenome}/{sample}/")
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['get_fastq_pe']['threads']
    benchmark:
        "benchmarks/{refGenome}/getfastq/{sample}_{run}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['get_fastq_pe']['mem'],
        tmpdir = get_big_temp
    shell:
        """
        set +e
        #delete existing prefetch file in case of previous run failure
        rm -rf {wildcards.run}
        ##attempt to get SRA file from NCBI (prefetch) or ENA (wget)
        prefetch --max-size 1T {wildcards.run}
        prefetchExit=$?
        if [[ $prefetchExit -ne 0 ]]
        then
            ffq --ftp {wildcards.run} | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | grep "fastq" | xargs curl --remote-name-all --output-dir {params.outdir}
        else
            fasterq-dump {wildcards.run} -O {params.outdir} -e {threads} -t {resources.tmpdir}
            pigz -p {threads} {params.outdir}{wildcards.run}*.fastq
        fi
        rm -rf {wildcards.run}
        """

rule fastp:
    input:
        unpack(get_reads)
    output:
        r1 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_1.fastq.gz",
        r2 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_2.fastq.gz",
        summ = "results/{refGenome}/summary_stats/{sample}/{run}.fastp.out"
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['fastp']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['fastp']['mem'],
    log:
        "logs/{refGenome}/fastp/{sample}/{run}.txt"
    benchmark:
        "benchmarks/{refGenome}/fastp/{sample}_{run}.txt"
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--thread {threads} "
        "--detect_adapter_for_pe "
        "-j /dev/null -h /dev/null "
        "2> {output.summ} > {log}"