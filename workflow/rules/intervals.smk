rule picard_intervals:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
        dictf = "results/{refGenome}/data/genome/{refGenome}.dict"
    output:
        intervals = temp("results/{refGenome}/intervals/picard_interval_list.list")
    params:
        minNmer = int(config['minNmer'])
    conda:
        '../envs/bam2vcf.yml'
    log:
        "logs/{refGenome}/picard_intervals/log.txt"
    benchmark:
        "benchmarks/{refGenome}/picard_intervals/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['process_ref']['mem']
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} OUTPUT_TYPE=ACGT &> {log}\n"

rule format_interval_list:
    input:
        intervals = "results/{refGenome}/intervals/picard_interval_list.list"
    output:
        intervals = "results/{refGenome}/intervals/master_interval_list.list"
    run:
        with open(output.intervals, "w") as out:
            with open(input.intervals, "r") as inp:
                for line in inp:
                    if not line.startswith("@"):
                        line = line.strip().split("\t")
                        chrom, start, end = line[0], line[1], line[2]
                        print(f"{chrom}:{start}-{end}", file=out)
    

checkpoint create_db_intervals:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        intervals = "results/{refGenome}/intervals/master_interval_list.list"
    output:
        out_dir = directory("results/{refGenome}/intervals/db_intervals"),
    params:
        max_intervals = get_db_interval_count
    log:
        "logs/{refGenome}/db_intervals/log.txt"
    benchmark:
        "benchmarks/{refGenome}/db_intervals/benchmark.txt"
    conda:
        '../envs/bam2vcf.yml'
    shell:
        """
        gatk SplitIntervals -L {input.intervals} \
        -O {output} -R {input.ref} -scatter {params} \
        -mode INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY &> {log}
        """

checkpoint create_gvcf_intervals:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        intervals = "results/{refGenome}/intervals/master_interval_list.list"
    output:
        out_dir = directory("results/{refGenome}/intervals/gvcf_intervals"),
    params:
        max_intervals = config["num_gvcf_intervals"]
    log:
        "logs/{refGenome}/gvcf_intervals/log.txt"
    benchmark:
        "benchmarks/{refGenome}/gvcf_intervals/benchmark.txt"
    conda:
        '../envs/bam2vcf.yml'
    shell:
        """
        gatk SplitIntervals -L {input.intervals} \
        -O {output} -R {input.ref} -scatter {params} \
        -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY  &> {log}
        """