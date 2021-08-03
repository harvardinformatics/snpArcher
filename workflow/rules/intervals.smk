rule process_ref:
    """
    This rule generates a .fai file from the reference genome, which are required for GATK to run properly. GATK also needs a .dict file, but this was previously generated.
    """
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output: 
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
    conda:
        "../envs/bam2vcf.yml"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "samtools faidx {input.ref} --output {output.fai}\n"
        "picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.dictf}\n"

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
        "../envs/bam2vcf.yml"
    log:
        "log/{Organism}/{refGenome}/picard_intervals/log"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} > {log}\n" 

rule create_intervals:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"
    params:
        maxIntervalLen = int(config['maxIntervalLen']),
        maxBpPerList = int(config['maxBpPerList']),
        maxIntervalsPerList = int(config['maxIntervalsPerList']),
        minNmer = int(config['minNmer'])
    output: 
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        LISTS = helperFun.createListsGetIndices(params.maxIntervalLen, params.maxBpPerList, params.maxIntervalsPerList, params.minNmer, config["output"], config["intDir"], wildcards, input.dictf, input.intervals)
