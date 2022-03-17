localrules: collect_fastp_stats, collect_sumstats
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
        alnSumMetsFiles = get_aln_sum_metrics,
        fastpFiles = get_fastpFiles,
        coverageFiles = get_coverageFiles
    output:
        config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_bam_sumstats.txt"
    run:
        FractionReadsPassFilter, NumReadsPassFilter = helperFun.collectFastpOutput(input.fastpFiles)
        aln_metrics = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
        SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
        helperFun.printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0])
