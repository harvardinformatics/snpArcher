rule genmap:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        bg = temp("results/{refGenome}/genmap/{refGenome}.genmap.bedgraph"),
        sorted_bg = "results/{refGenome}/genmap/sorted_mappability.bg"
    params:
        indir = os.path.join(workflow.default_remote_prefix, "results/{refGenome}/genmap_index"),
        outdir = os.path.join(workflow.default_remote_prefix, "results/{refGenome}/genmap")
    log:
        "logs/{refGenome}/genmap/log.txt"
    benchmark:
        "benchmarks/{refGenome}/genmap/benchmark.txt"
    conda:
        "../envs/genmap.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['genmap']['mem']
    threads:
        config['genmap']['threads'] 
    shell:
        # snakemake creates the output directory before the shell command, but genmap doesnt like this. so we remove the directory first.
        """
        rm -rf {params.indir} && genmap index -F {input.ref} -I {params.indir} &> {log}
        genmap map -K 150 -E 0 -I {params.indir} -O {params.outdir} -bg -T {threads} -v &> {log}
        sort -k1,1 -k2,2n {output.bg} > {output.sorted_bg} 2>> {log}
        """

rule genome_prep:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        twobit = "results/{refGenome}/callable_sites/genome.2bit",
        chrom = "results/{refGenome}/callable_sites/genome.sizes"
    conda:
        "../envs/callable.yml"
    log:
        "logs/{refGenome}/genome_prep/log.txt"
    benchmark:
        "benchmarks/{refGenome}/genome_prep/benchmark.txt"
    shell:
        "faToTwoBit {input.ref} {output.twobit} &> {log"
        "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}"

rule compute_d4:
    input:
        bam = "results/{refGenome}/bams/{sample}_final.bam",
        bai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    output:
        "results/{refGenome}/callable_sites/{sample}.mosdepth.global.dist.txt",
        temp("results/{refGenome}/callable_sites/{sample}.per-base.d4"),
        summary="results/{refGenome}/callable_sites/{sample}.mosdepth.summary.txt"
    conda:
        "../envs/callable.yml"
    log:
        "logs/{refGenome}/compute_d4/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/compute_d4/{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['compute_d4']['mem']
    threads:
        config['compute_d4']['threads']
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
        "../envs/callable.yml"
    log:
        "logs/{refGenome}/merge_d4/log.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_d4/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['merge_d4']['mem']
    shell:
        "d4tools merge {input.d4files} {output} &> {log}"

rule compute_covstats:
    input:
        d4 = "results/{refGenome}/callable_sites/all_samples.d4"
    output:
        stats = "results/{refGenome}/callable_sites/all_samples.covstats.txt"
    log:
        "logs/{refGenome}/compute_covstats/log.txt"
    benchmark:
        "benchmarks/{refGenome}/compute_covstats/benchmark.txt"
    conda:
        "../envs/callable.yml"
    shell:
        "d4tools stat {input.d4} | awk '{{ for(i=4; i<=NF;i++) j+=$i; print $1,$3,j; j=0 }}' > {output} 2> {log}"

rule create_cov_bed:
    input:
        stats = "results/{refGenome}/callable_sites/all_samples.covstats.txt",
        d4 = "results/{refGenome}/callable_sites/all_samples.d4"
    output:
        covbed = temp("results/{refGenome}/callable_sites/callable_sites_cov.bed")
    params:
        cov_threshold = config['cov_threshold']
    conda:
        "../envs/callable.yml"
    script:
        "../scripts/create_coverage_bed.py"

rule callable_bed:
    input:
        cov = "results/{refGenome}/callable_sites/callable_sites_cov.bed",
        map = "results/{refGenome}/genmap/sorted_mappability.bg"
    output:
        callable_sites = "result/{refGenome}/callable_sites.bed",
        tmp_cov = temp("results/{refGenome}/callable_sites/temp_cov.bed"),
        tmp_map = temp("results/{refGenome}/callable_sites/temp_map.bed")
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['callable_bed']['mem']
    params:
        merge = config['callable_merge'],
        mappability = config['mappability_min'],
    shell:
        """
        bedtools merge -i {input.cov} > {output.tmp_cov}
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.mappability}) print $1,$2,$3 }}' {input.map} > {output.tmp_map}
        bedtools intersect -a {output.tmp_cov} -b {output.tmp_map} | bedtools sort -i - | bedtools merge -d {params.merge} -i - > {output.callable_sites}
        """
