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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['compute_d4']['mem']
    threads:
        resources['compute_d4']['threads']
    params:
        prefix = os.path.join(workflow.default_remote_prefix, "results/{refGenome}/callable_sites/{sample}")
    shell:
        "mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}"

rule merge_d4:
    input:
        unpack(get_input_for_coverage)
    output:
        "results/{refGenome}/callable_sites/all_samples.d4"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{refGenome}/merge_d4/log.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_d4/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['merge_d4']['mem']
    shell:
        "d4tools merge {input.d4files} {output} &> {log}"

rule compute_covstats:
    input:
        d4 = "results/{refGenome}/callable_sites/all_samples.d4"
    output:
        stats = "results/{refGenome}/callable_sites/all_samples.covstats.txt",
        tmp_d4 = temp("results/{refGenome}/callable_sites/all_samples.covstats.tmp")
    log:
        "logs/{refGenome}/compute_covstats/log.txt"
    benchmark:
        "benchmarks/{refGenome}/compute_covstats/benchmark.txt"
    conda:
        "../envs/cov_filter.yml"
    shell:
        """
        d4tools stat {input.d4} > {output.tmp_d4} 2> {log}
        awk '{{ for(i=4; i<=NF;i++) j+=$i; print $1,$3,j; j=0 }}' {output.tmp_d4} > {output} 2>> {log}
        """

rule create_cov_bed:
    input:
        stats = "results/{refGenome}/callable_sites/all_samples.covstats.txt",
        d4 = "results/{refGenome}/callable_sites/all_samples.d4"
    output:
        covbed = temp("results/{refGenome}/callable_sites/callable_sites_cov.bed")
    params:
        cov_threshold = config['cov_threshold']
    conda:
        "../envs/cov_filter.yml"
    script:
        "../scripts/create_coverage_bed.py"

rule callable_bed:
    input:
        cov = "results/{refGenome}/callable_sites/callable_sites_cov.bed",
        map = "results/{refGenome}/callable_sites/callable_sites_map.bed"
    output:
        callable_sites = "results/{refGenome}/{prefix}_callable_sites.bed",
        tmp_cov = temp("results/{refGenome}/callable_sites/temp_cov.bed")
    conda:
        "../envs/cov_filter.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['callable_bed']['mem']
    params:
        merge = config['cov_merge']
    shell:
        """
        bedtools merge -i {input.cov} > {output.tmp_cov}
        bedtools sort -i {output.tmp_cov} | bedtools merge -d {params.merge} -i - > {output.tmp_cov}
        bedtools intersect -a {output.tmp_cov} -b {input.map} | bedtools sort -i - | bedtools merge -i - > {output.callable_sites}
        """
