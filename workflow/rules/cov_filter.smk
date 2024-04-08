rule compute_d4:
    input:
        bam = "results/{refGenome}/bams/{sample}_final.bam",
        bai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    output:
        "results/{refGenome}/callable_sites/{sample}.mosdepth.global.dist.txt",
        temp("results/{refGenome}/callable_sites/{sample}.per-base.d4"),
        summary="results/{refGenome}/callable_sites/{sample}.mosdepth.summary.txt"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{refGenome}/compute_d4/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/compute_d4/{sample}.txt"
    params:
        prefix = os.path.join(DEFAULT_STORAGE_PREFIX, "results/{refGenome}/callable_sites/{sample}")
    shell:
        "mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}"

rule merge_d4:
    input:
        ancient(unpack(get_input_for_coverage))
    output:
        "results/{refGenome}/callable_sites/all_samples.d4"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{refGenome}/merge_d4/log.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_d4/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['merge_d4']['mem'],
        machine_type = resources['sentieon_combine_gvcf']['machine_type'],
        disk_mb = resources['sentieon_combine_gvcf']['disk_mb']
    shell:
        "d4tools merge {input.d4files} {output} &> {log}"

rule collect_covstats:
    input:
        ancient(unpack(get_input_covstats))
    output:
        "results/{refGenome}/summary_stats/all_cov_sumstats.txt"  
    run:
        covStats = collectCovStats(input.covStatFiles)
        with open(output[0], "w") as f:
            print("chrom\tmean_cov\tstdev_cov", file=f)
            for chrom in covStats:
                print(chrom, covStats[chrom]['mean'], covStats[chrom]['stdev'], sep="\t", file=f)

rule create_cov_bed:
    input:
        stats = "results/{refGenome}/summary_stats/all_cov_sumstats.txt",
        d4 = "results/{refGenome}/callable_sites/all_samples.d4"
    output:
        covbed = "results/{refGenome}/callable_sites/{prefix}_callable_sites_cov.bed"
    benchmark:
        "benchmarks/{refGenome}/covbed/{prefix}_benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['callable_bed']['mem'],
        disk_mb = resources['callable_bed']['disk_mb']  
    params:
        cov_threshold_stdev = config["cov_threshold_stdev"],
        cov_threshold_lower = config["cov_threshold_lower"],
        cov_threshold_upper = config["cov_threshold_upper"],
        cov_threshold_rel = config["cov_threshold_rel"]
    conda:
        "../envs/cov_filter.yml"
    script:
        "../scripts/create_coverage_bed.py"

rule callable_bed:
    input:
        cov = "results/{refGenome}/callable_sites/{prefix}_callable_sites_cov.bed",
        map = "results/{refGenome}/callable_sites/{prefix}_callable_sites_map.bed"
    output:
        callable_sites = "results/{refGenome}/{prefix}_callable_sites.bed",
        tmp_cov = temp("results/{refGenome}/callable_sites/{prefix}_temp_cov.bed"),
    conda:
        "../envs/cov_filter.yml"
    benchmark:
        "benchmarks/{refGenome}/callable_bed/{prefix}_benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['callable_bed']['mem'],
        disk_mb = resources['callable_bed']['disk_mb']
    params:
        merge = config['cov_merge']
    shell:
        """
        bedtools sort -i {input.cov} | bedtools merge -d {params.merge} -i - > {output.tmp_cov}
        bedtools intersect -a {output.tmp_cov} -b {input.map} | bedtools sort -i - | bedtools merge -i - > {output.callable_sites}
        """
#awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' {output.callable_sites} > {output.bed_l}
#bed_l = temp("results/{refGenome}/callable_sites/{prefix}_length.txt")
