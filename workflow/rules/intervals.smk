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

checkpoint create_db_intervals:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "intervals.list"
    output:
        out_dir = directory(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "db_intervals"),
        
    params:
        max_intervals = get_db_interval_count
    conda:
        '../envs/bam2vcf.yml'
    shell:
        """
        gatk SplitIntervals -L {input.intervals} \
        -O {output} -R {input.ref} -scatter {params} \
        -mode INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY
        """

checkpoint create_gvcf_intervals:
    input:
        in_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"
    output:
        out_dir = directory(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "gvcf_intervals"),
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "intervals.list"
    params:
        max_intervals = config["maxNumIntervals"]
    script:
        "../scripts/make_intervals.py"