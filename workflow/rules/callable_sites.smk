import gzip
import json

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
      "faToTwoBit {input.ref} {output.twobit}\n"
      "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}"

rule bedgraphs:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix']
    output:
        temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}.sorted.bg")
    conda:
        "../envs/callable.yml"
    benchmark:
        "benchmarks/{Organism}/bedgraphs/{refGenome}_{Organism}_{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedtools']['mem']
    shell:
        "bedtools genomecov -ibam {input.bam} -bga | sort -k1,1 -k2,2n - > {output}"

rule merge_bedgraph:
    input:
        unpack(get_input_for_coverage)
    output:
        merge = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{Organism}.merge.bg")
    benchmark:
        "benchmarks/{Organism}/merge_begraph/{refGenome}_{Organism}.txt"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedtools']['mem']
    shell:
        "bedtools unionbedg -header -empty -g {input.chrom} -i {input.bedgraphs} > {output.merge}"

rule gzip_bedgraph:
    input:
        get_bedgraph_to_convert
    output:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz",
        idx = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gzi"
    conda:
        "../envs/callable.yml"
    shell:
        "bgzip -i -I {output.idx} -c {input} > {output.bgz}"

rule compute_covstats:
    input:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.bg.gz",
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.txt"
    run:
        covtot = 0.0
        meantot = 0.0
        bptot = 0.0
        with gzip.open(input.bgz, 'rt') as f:
            with gzip.open(output.cov, 'wt') as covbed:
                for line in f:
                    fields = line.split()
                    if (len(fields) < 4):
                        continue
                    try:
                        cov_fields = map(float, fields[3:])
                    except:
                        continue
                    covsum = sum(cov_fields)
                    mean = covsum / len(fields[3:])
                    bp = float(fields[2]) - float(fields[1])
                    covtot = covtot + (covsum * bp)
                    meantot = meantot + (mean * bp)
                    bptot = bptot + bp
                    print(fields[0], fields[1], fields[2], covsum, mean, file=covbed, sep="\t")
        stats = {"mean_summed_cov": covtot/bptot, "mean_ind_cov": meantot / bptot}
        with open(output.stats, 'w') as st:
            print(json.dumps(stats), file=st)

rule filter_bed:
    input:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.bg.gz",
        map = config['output'] + "{refGenome}/" + "genmap/{refGenome}.sorted_genmap.bg",
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.txt"
    output:
        callable_cov = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_cov.bed"),
        callable_map = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_map.bed")
    params:
        mappability = config['mappability_min'],
        low_cov = config['low_cov_thresh'],
        high_cov = config['high_cov_thresh'],
        cov_type = config['cov_type']
    run:
        stats = {}
        with open(input.stats) as stats:
            stats = json.load(stats)
        min_cov = float(stats[params.cov_type]) * params.low_cov
        max_cov = float(stats[params.cov_type]) * params.high_cov
        record_num = 3 if params.cov_type == "mean_summed_cov" else 4
        with gzip.open(input.cov, 'rt') as cov:
            with open(output.callable_cov, 'w') as cov_out:
                for line in cov:
                    fields = line.split()
                    if float(fields[record_num]) >= min_cov and float(fields[record_num]) <= max_cov:
                        print(fields[0], fields[1], fields[2], sep="\t", file=cov_out)
        with open(input.map) as map:
            with open(output.callable_map, 'w') as map_out:
                for line in map:
                    fields=line.split()
                    if float(fields[3]) >= params.mappability:
                        print(fields[0], fields[1], fields[2], sep="\t", file=map_out)

rule callable_bed:
    input:
        callable_cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_cov.bed",
        callable_map = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites_map.bed"
    output:
        callable_sites = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".callable_sites.bed"
    conda:
        "../envs/callable.yml"
    params:
        merge = config['callable_merge']
    shell:
        "bedtools intersect -a {input.callable_cov} -b {input.callable_map} | bedtools merge -d {params.merge} -i - > {output.callable_sites}"
