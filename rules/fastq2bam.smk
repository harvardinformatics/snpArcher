localrules: collect_sumstats


rule index_ref:
    input:
        ref = config['ref']
    output: 
       config['ref'] + ".sa",
       config['ref'] + ".pac",
       config['ref'] + ".bwt",
       config['ref'] + ".ann",
       config['ref'] + ".amb"
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['index_ref']['mem'] 
    shell:
        "bwa index {input.ref}"

rule fastp:
    input:
        r1 = fastqDir + "{sample}" + fastq_suffix1,
        r2 = fastqDir + "{sample}" + fastq_suffix2
    output: 
        r1 = fastqFilterDir + "{sample}_fastp" + fastq_suffix1,
        r2 = fastqFilterDir + "{sample}_fastp" + fastq_suffix2,
        summ = sumstatDir + "{sample}_fastp.out"
    conda:
        "../envs/fastq2bam.yml"
    threads: res_config['fastp']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['fastp']['mem'] 
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--thread {threads} "
        "--detect_adapter_for_pe "
        "2> {output.summ}"

rule bwa_map:
    input:
        ref = config['ref'],
        r1 = fastqFilterDir + "{sample}_fastp" + fastq_suffix1,
        r2 = fastqFilterDir + "{sample}_fastp" + fastq_suffix2,
        # the following files are bwa index files that aren't directly input into command below, but needed
        sa = config['ref'] + ".sa",
        pac = config['ref'] + ".pac",
        bwt = config['ref'] + ".bwt",
        ann = config['ref'] + ".ann",
        amb = config['ref'] + ".amb"
    output: 
        bamDir + "{sample}.bam"
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    conda:
        "../envs/fastq2bam.yml"
    threads: res_config['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bwa_map']['mem'] 
    shell:
        "bwa mem -M -t {threads} -R \'{params.rg}\' {input.ref} {input.r1} {input.r2} | "
        "samtools view -Sb - > {output}"

rule sort_bam:
    input: 
        bamDir + "{sample}.bam"
    output: 
        temp(bamDir + "{sample}_sorted.bam"),
        temp(bamDir + "{sample}_sorted.bai")
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['sort_bam']['mem'] 
    shell:
        "picard SortSam I={input} O={output[0]} SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR={bamDir}"

rule dedup:
    input: 
        bamDir + "{sample}_sorted.bam",
        bamDir + "{sample}_sorted.bai"
    output:
        dedupBam = bamDir + "{sample}_dedup.bam",
        dedupMet = sumstatDir + "{sample}_dedupMetrics.txt",
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['dedup']['mem'] 
    shell:
        "picard MarkDuplicates I={input[0]} O={output.dedupBam} METRICS_FILE={output.dedupMet} REMOVE_DUPLICATES=false TAGGING_POLICY=All\n"
        "picard BuildBamIndex I={output.dedupBam} "

rule bam_sumstats:
    input: 
        bam = bamDir + "{sample}_dedup.bam",
        ref = config['ref']

    output: 
        cov = sumstatDir + "{sample}_coverage.txt",
        alnSum = sumstatDir + "{sample}_AlnSumMets.txt",
        val = sumstatDir + "{sample}_validate.txt"
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam_sumstats']['mem'] 
    shell:
        "samtools coverage --output {output.cov} {input.bam}\n"
        "picard CollectAlignmentSummaryMetrics I={input.bam} R={input.ref} O={output.alnSum}\n"
        # The following ValidateSamFile exits with non-zero status when a BAM file contains errors, 
        # causing snakemake to exit and remove these output files.  I cirumvent this by appending "|| true".
        # I also ignore "INVALID_TAG_NM" because it isn't used by GATK but causes errors at this step
        "picard ValidateSamFile I={input.bam} R={input.ref} O={output.val} IGNORE=INVALID_TAG_NM || true"
		
rule collect_sumstats:
    input:
        fastpFiles = expand(sumstatDir + "{sample}_fastp.out", sample=SAMPLES),
        dedupFiles = expand(sumstatDir + "{sample}_dedupMetrics.txt", sample=SAMPLES),
        alnSumMetsFiles = expand(sumstatDir + "{sample}_AlnSumMets.txt", sample=SAMPLES),
        coverageFiles = expand(sumstatDir + "{sample}_coverage.txt", sample=SAMPLES),
        validateFiles = expand(sumstatDir + "{sample}_validate.txt", sample=SAMPLES)
    output:
         config["fastq2bamDir"] + "bam_sumstats.txt"
    run:
        FractionReadsPassFilter, NumFilteredReads = helperFun.collectFastpOutput(input.fastpFiles)
        PercentDuplicates = helperFun.collectDedupMetrics(input.dedupFiles)
        PercentHQreads, PercentHQbases = helperFun.collectAlnSumMets(input.alnSumMetsFiles)
        SeqDepths, CoveredBases = helperFun.collectCoverageMetrics(input.coverageFiles)
        validateSams = helperFun.collectValidationStatus(input.validateFiles)

        helperFun.printBamSumStats(FractionReadsPassFilter, NumFilteredReads, PercentDuplicates, PercentHQreads, PercentHQbases, SeqDepths, CoveredBases, validateSams, config["fastq2bamDir"])

