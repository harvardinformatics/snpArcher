rule genmap_index:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
    log:
        "logs/{Organism}/{refGenome}/genmap_index/log"
    conda:
        "../envs/genmap.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
    output:
        temp(directory(config['output'] + "{Organism}/{refGenome}/" + "genmap/genmap.index"))
    shell:
        "genmap index -F {input.ref} -I {output} &> {log}"

rule genmap_map:
    input:
        config['output'] + "{Organism}/{refGenome}/" + "genmap/genmap.index"
    log:
        "logs/{Organism}/{refGenome}/genmap_map/log"
    params:
        outdir = config['output'] + "{Organism}/{refGenome}/" + "genmap"
    conda:
        "../envs/genmap.yml"
    threads:
        res_config['genmap']['threads']
    output:
        bg = temp(config['output'] + "{Organism}/{refGenome}/" + "genmap/genmap.bedgraph")
    shell:
        "genmap map -K 150 -E 0 -I {input} -O {params.outdir} -bg -T {threads} -v  > {log}"

rule sort_genmap:
    input:
        config['output'] + "{Organism}/{refGenome}/" + "genmap/genmap.bedgraph"
    output:
        config['output'] + "{Organism}/{refGenome}/" + "genmap/sorted_genmap.bg"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule picard_intervals:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
    output:
        intervals = temp(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list")
    params:
        minNmer = int(config['minNmer'])
    conda:
        '../envs/bam2vcf.yml'
    log:
        "logs/{Organism}/{refGenome}/picard_intervals/log"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} &> {log}\n"

checkpoint create_intervals:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"
    params:
        maxIntervalLen = int(config['maxIntervalLen']),
        maxBpPerList = int(config['maxBpPerList']),
        maxIntervalsPerList = int(config['maxIntervalsPerList']),
        minNmer = int(config['minNmer']),
        max_intervals = config['maxNumIntervals']
    output:
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem']
    run:
        if config['split_by_n']:
            LISTS = helperFun.createListsGetIndices(params.maxIntervalLen, params.maxBpPerList, params.maxIntervalsPerList, params.minNmer, config["output"], config["intDir"], wildcards, input.dictf, input.intervals)
        else:
            LISTS = make_intervals(config["output"], config["intDir"], wildcards, input.dictf, params.max_intervals)
