import gzip
import json

## RULES ##

## plan is to make a callable sites that filters on both mappability and coverage


rule compute_covstats:
    input:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.bg.gz",
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "covstats.txt"
    run:
        covtot = 0
        meantot = 0
        bptot = 0
        with gzip.open(input.bgz, 'rt') as f:
            with gzip.open(output.cov, 'wt') as covbed:
                for line in f:
                    fields = line.split()
                    if (len(fields) < 4):
                        continue
                    try:
                        cov_fields = map(int, fields[3:])
                    except:
                        continue
                    covsum = sum(cov_fields)
                    mean = covsum / len(fields[3:])
                    bp = int(fields[2]) - int(fields[1])
                    covtot = covtot + (covsum * bp)
                    meantot = meantot + (mean * bp)
                    bptot = bptot + bp
                    print(fields[0], fields[1], fields[2], covsum, mean, file=covbed, sep="\t")
        stats = {"mean_summed_cov": covtot/bptot, "mean_ind_cov": meantot / bptot}
        with open(output.stats) as st:
            print(json.dumps(stats), file=st)

rule filter_bed:
    input:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.bg.gz",
        map = config['output'] + "{refGenome}/" + "genmap/{refGenome}.sorted_genmap.bg",
        stats = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "covstats.txt"
    output:
        callable_cov = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "callable_sites_cov.bed"),
        callable_map = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "callable_sites_map.bed")
    params:
        mappability = 1,
        low_cov = 0.5,
        high_cov = 2
    run:
        stats = {}
        with open(input.stats) as stats:
            stats = json.load(stats)
        min_cov = stats['mean_summed_cov'] * params.low_cov
        max_cov = stats['mean_summed_cov'] * params.high_cov
        with gzip.open(input.cov, 'rt') as cov:
            with open(output.callable_cov) as cov_out:
                for line in cov:
                    fields = line.split()
                    if fields[3] >= min_cov and fields[3] <= max_cov:
                        print(fields[0], fields[1], fields[2], sep="\t", file=cov_out)
        with open(input.map) as map:
            with open(output.callable_map) as map_out:
                for line in map:
                    fields=line.split()
                    if fields[3] >= mappability:
                        print(fields[0], fields[1], fields[2], sep="\t", file=map_out)

rule callable_bed:
    input:
        callable_cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "callable_sites_cov.bed",
        callable_map = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "callable_sites_map.bed"
    output:
        callable_sites = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "callable_sites.bed"
    envs:
        "../envs/callable.yml"
    params:
        merge = 100
    shell:
        bedtools intersect -a {input.callable_cov} -b {input.callable_map} | bedtools merge -d {params.merge} -i - > {output.callable_sites}
