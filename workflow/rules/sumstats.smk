rule bam_sumstats:
    input:
        unpack(get_bams),
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        cov = "results/{refGenome}/summary_stats/{sample}_coverage.txt",
        alnSum = "results/{refGenome}/summary_stats/{sample}_AlnSumMets.txt",
    conda:
        "../envs/fastq2bam.yml"
    shell:
        """
        samtools coverage --output {output.cov} {input.bam}
        samtools flagstat -O tsv {input.bam} > {output.alnSum}
        """

rule sentieon_bam_stats:
    input:
        unpack(get_bams),
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        ref = "results/{refGenome}/data/genome/{refGenome}.fna"
    params:
        lic = config['sentieon_lic']
    output:
        insert_file = "results/{refGenome}/summary_stats/{sample}_insert_metrics.txt",
        qd = "results/{refGenome}/summary_stats/{sample}_qd_metrics.txt",
        gc = "results/{refGenome}/summary_stats/{sample}_gc_metrics.txt",
        gc_summary = "results/{refGenome}/summary_stats/{sample}_gc_summary.txt",
        mq = "results/{refGenome}/summary_stats/{sample}_mq_metrics.txt"
    conda:
        "../envs/sentieon.yml"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} \
        -t {threads} -i {input.bam} \
        --algo MeanQualityByCycle {output.mq} \
        --algo QualDistribution {output.qd} \
        --algo GCBias --summary {output.gc_summary} {output.gc} \
        --algo InsertSizeMetricAlgo {output.insert_file}
        """

rule collect_fastp_stats:
    input:
        collect_fastp_stats_input
    output:
        "results/{refGenome}/summary_stats/{sample}_fastp.out"
    run:
        combine_fastp_files(input, output)

rule collect_sumstats:
    input:
        unpack(get_input_sumstats)
    output:
        "results/{refGenome}/summary_stats/{prefix}_bam_sumstats.txt"
    run:
        if not config['sentieon']:
            FractionReadsPassFilter, NumReadsPassFilter = collectFastpOutput(input.fastpFiles)
            aln_metrics = collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = collectCoverageMetrics(input.coverageFiles)
            printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0])
        else:
            FractionReadsPassFilter, NumReadsPassFilter = collectFastpOutput(input.fastpFiles)
            aln_metrics = collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = collectCoverageMetrics(input.coverageFiles)
            median_inserts, median_insert_std = collect_inserts(input.insert_files)
            printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0], median_inserts, median_insert_std)