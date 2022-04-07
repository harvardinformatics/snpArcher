localrules: genome_prep

## RULES ##

rule genmap_index:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
    log:
        "logs/{refGenome}/genmap_index.log"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
    output:
        temp(directory(config['output'] + "{refGenome}/" + "genmap.index"))
    shell:
        "genmap index -F {input.ref} -I {output} &> {log}"

rule genmap_map:
    input:
        config['output'] + "{refGenome}/" + "genmap.index"
    log:
        "logs/{refGenome}/genmap_map.log"
    params:
        outdir = config['output'] + "{refGenome}/" + "genmap"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
    threads:
        res_config['genmap']['threads']
    output:
        temp(config['output'] + "{refGenome}/" + "genmap/{refGenome}.genmap.bedgraph")
    shell:
        "genmap map -K 150 -E 0 -I {input} -O {params.outdir} -bg -T {threads} -v  > {log}"

rule sort_genmap:
    input:
        config['output'] + "{refGenome}/" + "genmap/{refGenome}.genmap.bedgraph"
    output:
        config['output'] + "{refGenome}/" + "genmap/{refGenome}.sorted_genmap.bg"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap_sort']['mem']
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule genome_prep:
  input:
      ref = config["refGenomeDir"] + "{refGenome}.fna"
  output:
      twobit = config['output'] + "{refGenome}/" + "{refGenome}" + ".2bit",
      chrom = config['output'] + "{refGenome}/" + "{refGenome}" + ".sizes"
  conda:
      "../envs/callable.yml"
  shell:
      "faToTwoBit {input.ref} {output.twobit}"
      "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}"

rule compute_d4:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix']
    output:
        config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}.mosdepth.global.dist.txt",
        temp(config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}.per-base.d4"),
        summary=config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}.mosdepth.summary.txt"
    conda:
        "../envs/callable.yml"
    benchmark:
        "benchmarks/{Organism}/bedgraphs/{refGenome}_{Organism}_{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['compute_d4']['mem']
    threads:
        res_config['compute_d4']['threads']
    params:
        prefix = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}"
    shell:
        "mosdepth --d4 -t {threads} {params.prefix} {input.bam}"

rule merge_d4:
    input:
        unpack(get_input_for_coverage)
    output:
        config['output'] + "{Organism}/{refGenome}/{refGenome}_{Organism}.d4"
    benchmark:
        "benchmarks/{Organism}/merge_d4/{refGenome}_{Organism}.txt"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['merge_d4']['mem']
    shell:
        "d4tools merge {input.d4files} {output}"

rule compute_covstats:
    input:
        d4 = config['output'] + "{Organism}/{refGenome}/{refGenome}_{Organism}.d4"
    output:
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.txt"
    conda:
        "../envs/callable.yml"
    shell:
        "d4tools stat {input.d4} | awk '{{ for(i=4; i<=NF;i++) j+=$i; print $1,$3,j; j=0 }}' > {output}"

rule create_cov_bed:
    input:
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.txt",
        d4 = config['output'] + "{Organism}/{refGenome}/{refGenome}_{Organism}.d4"
    output:
        covbed = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_cov.bed")
    params:
        cov_threshold = config['cov_threshold']
    conda:
        "../envs/callable.yml"
    script:
        "../scripts/create_coverage_bed.py"

rule callable_bed:
    input:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_cov.bed",
        map = config['output'] + "{refGenome}/" + "genmap/{refGenome}.sorted_genmap.bg"
    output:
        callable_sites = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites.bed"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['callable_bed']['mem']
    params:
        merge = config['callable_merge'],
        mappability = config['mappability_min'],
        #possibly better ways to create temp files in snakemake but this will work for now
        TEMP_cov = config['tmp_dir'] + "{Organism}_{refGenome}_TEMP.cov.bed",
        TEMP_map = config['tmp_dir'] + "{Organism}_{refGenome}_TEMP.map.bed"
    shell:
        "bedtools merge {input.callable_cov} > {params.TEMP_cov}"
        """awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.mappability}) print $1,$2,$3 }}' {input.map} > {params.TEMP_map}"""
        "bedtools intersect -a {params.TEMP_cov} -b {params.TEMP_map} | bedtools sort -i - | bedtools merge -d {params.merge} -i - > {output.callable_sites}"
        "rm {params.TEMP_cov}"
        "rm {params.TEMP_map}"
