localrules: collect_sumstats
rule bam_sumstats:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam",
        bai = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam.bai",
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
rule sentieon_bam_stats:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam",
        bai = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam.bai",
        indices = ancient(expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["fai", "sa", "pac", "bwt", "ann", "amb"])),
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        lic = ancient(config['sentieon_lic'])
    output:
        insert_file = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_insert_metrics.txt",
        qd = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_qd_metrics.txt",
        gc = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_gc_metrics.txt",
        gc_summary = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_gc_summary.txt",
        mq = config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_mq_metrics.txt"
    conda:
        "../envs/sentieon.yml"
    threads: 8
    shell:
        """
        export SENTIEON_LICENSE={input.lic}
        sentieon driver -r {input.ref} \
        -t {threads} -i {input.bam} \
        --algo MeanQualityByCycle {output.mq} \
        --algo QualDistribution {output.qd} \
        --algo GCBias --summary {output.gc_summary} {output.gc} \
        --algo InsertSizeMetricAlgo {output.insert_file}
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
        unpack(get_input_sumstats)
    output:
        config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_bam_sumstats.txt"
    run:
        if not config['sentieon']:
            FractionReadsPassFilter, NumReadsPassFilter = helperFun.collectFastpOutput(input.fastpFiles)
            aln_metrics = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
            helperFun.printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0])
        else:
            FractionReadsPassFilter, NumReadsPassFilter = helperFun.collectFastpOutput(input.fastpFiles)
            aln_metrics = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
            median_inserts, median_insert_std = helperFun.collect_inserts(input.insert_files)
            helperFun.printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0], median_inserts, median_insert_std)