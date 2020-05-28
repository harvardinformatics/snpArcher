rule index_ref:
    input:
        ref = config['ref']
    output: 
       config['ref'] + ".sa",
       config['ref'] + ".pac",
       config['ref'] + ".bwt",
       config['ref'] + ".ann",
       config['ref'] + ".amb"
    shell:
        "bwa index {input.ref}"

rule bwa_map:
    input:
        ref = config['ref'],
        r1 = fastqDir + "{sample}" + fastq_suffix1,
        r2 = fastqDir + "{sample}" + fastq_suffix2,
        # the following files are bwa index files that aren't directly input into command below, but needed
        sa = config['ref'] + ".sa",
        pac = config['ref'] + ".pac",
        bwt = config['ref'] + ".bwt",
        ann = config['ref'] + ".ann",
        amb = config['ref'] + ".amb"
    output: 
        bamDir + "{sample}.bam"
    threads: 
        CLUSTER["bwa_map"]["n"]
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    shell:
        "bwa mem -M -t {threads} -R \'{params.rg}\' {input.ref} {input.r1} {input.r2} | "
        "/n/home11/bjarnold/programs/samtools-1.10/samtools view -Sb - > {output}"

rule sort_bam:
    input: 
        bamDir + "{sample}.bam"
    output: 
        temp(bamDir + "{sample}_sorted.bam"),
        temp(bamDir + "{sample}_sorted.bai")
    run:
        shell("java -jar /n/home11/bjarnold/picard.jar SortSam TMP_DIR={bamDir}tmp I={input} O={output[0]} SORT_ORDER=coordinate CREATE_INDEX=true")
        shell("rm {input}")

rule dedup:
    input: 
        bamDir + "{sample}_sorted.bam",
        bamDir + "{sample}_sorted.bai"
    output:
        dedupBam = temp(bamDir + "{sample}_dedup.bam"),
        dedupMet = sumstatDir + "{sample}_dedupMetrics.txt",
        dedupBamSort = bamDir + "{sample}_dedupSort.bam"
    run:
        shell("java -jar /n/home11/bjarnold/picard.jar MarkDuplicates TMP_DIR={bamDir}tmp I={input[0]} O={output.dedupBam} METRICS_FILE={output.dedupMet} REMOVE_DUPLICATES=false TAGGING_POLICY=All")
        shell("java -jar /n/home11/bjarnold/picard.jar SortSam TMP_DIR={bamDir}tmp I={output.dedupBam} O={output.dedupBamSort} SORT_ORDER=coordinate CREATE_INDEX=true")

rule bam_sumstats:
    input: 
        bam = bamDir + "{sample}_dedupSort.bam",
        ref = config['ref']

    output: 
        cov = sumstatDir + "{sample}_coverage.txt",
        alnSum = sumstatDir + "{sample}_AlnSumMets.txt",
        val = sumstatDir + "{sample}_validate.txt"
    run:
        shell("/n/home11/bjarnold/programs/samtools-1.10/samtools coverage --output {output.cov} {input.bam}")
        shell("java -jar /n/home11/bjarnold/picard.jar CollectAlignmentSummaryMetrics I={input.bam} R={input.ref} O={output.alnSum}")
        # The following ValidateSamFile exits with non-zero status when a BAM file contains errors, 
        # causing snakemake to exit and remove these output files.  I cirumvent this by appending "|| true".
        shell("java -jar /n/home11/bjarnold/picard.jar ValidateSamFile I={input.bam} R={input.ref} O={output.val} || true")
		
rule collect_sumstats:
    input:
        dedupFiles = expand(sumstatDir + "{sample}_dedupMetrics.txt", sample=SAMPLES),
        alnSumMetsFiles = expand(sumstatDir + "{sample}_AlnSumMets.txt", sample=SAMPLES),
        coverageFiles = expand(sumstatDir + "{sample}_coverage.txt", sample=SAMPLES),
        validateFiles = expand(sumstatDir + "{sample}_validate.txt", sample=SAMPLES)
    output:
        "bam_sumstats.txt"
    run:
        PercentDuplicates = helperFun.collectDedupMetrics(input.dedupFiles)
        PercentHQreads, PercentHQbases = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
        SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
        validateSams = helperFun.collectValidationStatus(input.validateFiles)

        helperFun.printBamSumStats(PercentDuplicates, PercentHQreads, PercentHQbases, SeqDepths, CoveredBases, validateSams)

