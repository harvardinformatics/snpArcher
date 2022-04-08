ruleorder: index_ref > download_reference
localrules: download_reference, index_ref

rule genmap_index:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
    log:
        "logs/{refGenome}/genmap_index.log"
    conda:
        "../envs/genmap.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
    params:
        os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap_index"))
    output:
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.{ext}.{ext2}", ext=['ids', 'info', 'txt'], ext2=['concat', 'limits']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext}.sbl", ext=['drp', 'drv']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext}.sbl", ext=['drp', 'drv']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.sa.{ext}", ext=['ind', 'len', 'val'])
    shell:
        # snakemake creates the output directory before the shell command, but genmap doesnt like this. so we remove the directory first.
        "rm -rf {params} && genmap index -F {input.ref} -I {params} &> {log}"

rule genmap_map:
    input:
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.{ext}.{ext2}", ext=['ids', 'info', 'txt'], ext2=['concat', 'limits']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext}.sbl", ext=['drp', 'drv']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext}.sbl", ext=['drp', 'drv']),
        expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.sa.{ext}", ext=['ind', 'len', 'val'])
    log:
        "logs/{refGenome}/genmap_map.log"
    params:
        indir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap_index")),
        outdir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap"))
    conda:
        "../envs/genmap.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
    threads:
        res_config['genmap']['threads']
    output:
        temp(config['output'] + "{refGenome}/" + "genmap/{refGenome}.genmap.bedgraph")
    shell:
        "genmap map -K 150 -E 0 -I {params.indir} -O {params.outdir} -bg -T {threads} -v  > {log}"

rule sort_genmap:
    input:
        config['output'] + "{refGenome}/" + "genmap/{refGenome}.genmap.bedgraph"
    output:
        config['output'] + "{refGenome}/" + "genmap/{refGenome}.sorted_genmap.bg"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['genmap_sort']['mem']
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
        "logs/{refGenome}/{Organism}.picard_intervals.log"
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
        max_intervals = config['maxNumIntervals'],
        missingBpTolerance = config['missingBpTolerance']
    output:
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem']
    run:
        if config['split_by_n']:
            LISTS = helperFun.createListsGetIndices(params.missingBpTolerance, params.maxIntervalLen, params.maxBpPerList, params.maxIntervalsPerList, params.minNmer, config["output"], config["intDir"], wildcards, input.dictf, input.intervals)
        else:
            LISTS = make_intervals(config["output"], config["intDir"], wildcards, input.dictf, params.max_intervals)
