localrules: collect_sumstats, download_reference
ruleorder: index_ref > download_reference

### RULES ###

rule get_fastq_pe:
    output:
        temp(config["fastqDir"] + "{Organism}/{sample}/{run}_1.fastq.gz"),
        temp(config["fastqDir"] + "{Organism}/{sample}/{run}_2.fastq.gz")
    params:
        outdir = config["fastqDir"] + "{Organism}/{sample}/",
        tmpdir = config['tmp_dir'],
        sra_url = lambda wildcards: get_ena_url(wildcards)["sra_url"],
        fastq_url = lambda wildcards: get_ena_url(wildcards)["fastq_url"]
    conda:
        "../envs/fastq2bam.yml"
    threads:
        res_config['get_fastq_pe']['threads']
    log:
        "logs/{Organism}/get_fastq/{sample}/{run}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['get_fastq_pe']['mem']
    shell:
        """
        set +e

        ##attempt to get SRA file from NCIB (prefetch) or ENA (wget)
        prefetch {wildcards.run}
        prefetchExit=$?
        if [[ $prefetchExit -ne 0 ]]
        then
            wget -O {wildcards.run} {params.sra_url}
        fi
        ##if this succeeded, we'll have the correct file in our working directory
        if [[ -s {wildcards.run} ]]
        then
            fasterq-dump {wildcards.run} -O {params.outdir} -e {threads} -t {params.tmpdir}
            pigz -p {threads} {params.outdir}{wildcards.run}*.fastq
        else
            wget -P {params.outdir} {params.fastq_url}/{wildcards.run}_1.fastq.gz
            wget -P {params.outdir} {params.fastq_url}/{wildcards.run}_2.fastq.gz
        fi
        rm -rf {wildcards.run}
        """

rule download_reference:
    output:
        outdir = directory(config["refGenomeDir"] + "{refGenome}"),
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    params:
        dataset = config["refGenomeDir"] + "{refGenome}_dataset.zip"
    log:
        "logs/dl_reference/{refGenome}.log"
    conda:
        "../envs/fastq2bam.yml"
    shell:
        "datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} &> {log}"
        "&& 7z x {params.dataset} -aoa -o{output.outdir}"
        "&& cat {output.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}"

rule index_ref:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        indexes = expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"]),
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict"
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['index_ref']['mem']
    log:
        "logs/index_ref/{refGenome}.log"
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} --output {output.fai}
        samtools dict {input.ref} -o {output.dictf} >> {log} 2>&1
        """

rule fastp:
    input:
        unpack(get_reads)
    output:
        r1 = temp(config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_1.fastq.gz"),
        r2 = temp(config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_2.fastq.gz"),
        summ = temp(config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}/{run}.out")
    conda:
        "../envs/fastq2bam.yml"
    threads:
        res_config['fastp']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['fastp']['mem']
    log:
        "logs/{Organism}/fastp/{refGenome}_{sample}_{run}.txt"
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--thread {threads} "
        "--detect_adapter_for_pe "
        "-j /dev/null -h /dev/null "
        "2> {output.summ} > {log}"

rule bwa_map:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        r1 = config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_1.fastq.gz",
        r2 = config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_2.fastq.gz",
        indices = expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"])
    output:
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}/{run}.bam")
    params:
        get_read_group
    conda:
        "../envs/fastq2bam.yml"
    threads:
        res_config['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bwa_map']['mem']
    log:
        "logs/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    benchmark:
        "benchmarks/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {threads} {params} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam} -"

rule merge_bams:
    input:
        lambda wildcards:
        expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    output:
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam"),
        bai = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['merge_bams']['mem']
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam}"

rule dedup:
    input:
        get_bams_for_dedup
    output:
        dedupBam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix'],
        dedupBai = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam.bai",
    conda:
        "../envs/sambamba.yml"
    resources:
        threads = res_config['dedup']['threads'],
        mem_mb = lambda wildcards, attempt: attempt * res_config['dedup']['mem']
    log:
        "logs/{Organism}/dedup/{refGenome}_{sample}.txt"
    benchmark:
        "benchmarks/{Organism}/dedup/{refGenome}_{sample}.txt"
    shell:
        "sambamba markdup -t {threads} {input} {output.dedupBam} 2> {log}"

rule bam_sumstats:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix'],
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_coverage.txt",
        alnSum = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_AlnSumMets.txt",
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam_sumstats']['mem']
    shell:
        """
        samtools coverage --output {output.cov} {input.bam}
        samtools flagstat -O tsv {input.bam} > {output.alnSum}
        """

rule collect_fastp_stats:
    input:
        lambda wildcards:
            expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{{sample}}/{run}.out", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    output:
        config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_fastp.out"
    shell:
        "cat {input} > {output}"

rule collect_sumstats:
    input:
        unpack(get_sumstats)
    output:
        config['output'] + "{Organism}/{refGenome}/" + "bam_sumstats.tsv"
    run:
        FractionReadsPassFilter, NumFilteredReads = helperFun.collectFastpOutput(input.fastpFiles)
        aln_metrics = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
        SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
        helperFun.printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumFilteredReads, output[0])
