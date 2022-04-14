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
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} OUTPUT_TYPE=ACGT &> {log}\n"

rule prepare_db_intervals:
    input:
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"
    output:
        intervals = temp((config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "db.interval_list"))
    params:
        max_interval_size = config["maxIntervalLen"]
    conda:
        '../envs/bam2vcf.yml'
    shell:
        "picard IntervalListTools I={input} BRK={params.max_interval_size} O={output}"

rule create_db_intervals:
    input:
        in_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "db.interval_list"
    output:
        out_files = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config["intDir"] + "db_interval_{i}.list", i=range(config["DBmaxNumIntervals"]))
    params:
        max_intervals = config["DBmaxNumIntervals"]
    script:
        "../scripts/make_intervals.py"

rule create_gvcf_intervals:
    input:
        in_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"
    output:
        out_files = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config["intDir"] + "gvcf_interval_{i}.list", i=range(config["maxNumIntervals"]))
    params:
        max_intervals = config["maxNumIntervals"]
    script:
        "../scripts/make_intervals.py"