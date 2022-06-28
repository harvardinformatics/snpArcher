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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['get_fastq_pe']['mem']
    shell:
        """
        set +e

        #delete existing prefetch file in case of previous run failure
        rm -rf {wildcards.run}

        ##attempt to get SRA file from NCBI (prefetch) or ENA (wget)
        prefetch --max-size 1T {wildcards.run}
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
    input:
        ref = get_ref
    output:
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    params:
        dataset = config["refGenomeDir"] + "{refGenome}_dataset.zip",
        outdir = config["refGenomeDir"] + "{refGenome}"
    conda:
        "../envs/fastq2bam.yml"
    shell:
        """
        if [ -z "{input.ref}" ]  # check if this is empty
        then
            datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} \
            && 7z x {params.dataset} -aoa -o{params.outdir} \
            && cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}
        else
            cp {input.ref} {output.ref}
        fi
        """
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
    benchmark:
        "benchmarks/{Organism}/fastp/{refGenome}_{sample}_{run}.txt"
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
        RG = get_read_group,
        bwa_threads = res_config['bwa_map']['threads'],
        sort_threads = res_config['sort_bam']['threads'],
        sort_mem = res_config['sort_bam']['mem_per_thread']
    conda:
        "../envs/fastq2bam.yml"
    threads:
        lambda wildcards: res_config['sort_bam']['threads'] + res_config['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * (res_config['bwa_map']['mem'] + (res_config['sort_bam']['mem_per_thread'] * res_config['sort_bam']['threads']))
    log:
        "logs/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    benchmark:
        "benchmarks/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {params.bwa_threads} {params.RG} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort --threads {params.sort_threads} -m {params.sort_mem}M -u - > {output.bam}"

rule merge_bams:
    input:
        lambda wildcards:
        expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    output:
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam"),
        bai = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    threads:
        res_config['merge_bams']['threads']
    params:
        comp_threads = lambda wildcards, threads: threads - 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['merge_bams']['mem']
    shell:
        "samtools merge -f --threads {params.comp_threads} -o {output.bam} {input} && samtools index -@ {params.comp_threads} {output.bam}"

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

rule bam2gvcf:
    """
    This rule scatters analyses over two dimensions: sample name and list file. For each BAM file, one per sample,
    a GVCF is created for all the scaffolds present in a given list file.
    """
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix'],
        l = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "gvcf_intervals/{list}.list"
    output:
        gvcf = temp(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "{list}.raw.g.vcf.gz"),
        gvcf_idx = temp(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "{list}.raw.g.vcf.gz.tbi"),
        
    resources:
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['bam2gvcf']['mem'] - 1000)  # this is the maximum amount given to java
    log:
        "logs/{Organism}/bam2gvcf/{refGenome}_{sample}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/bam2gvcf/{refGenome}_{sample}_{list}.txt"
    params:
        minPrun = config['minP'],
        minDang = config['minD'],
        
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk HaplotypeCaller "
        "--java-options \"-Xmx{resources.reduced}m\" "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-L {input.l} "
        "--emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang} &> {log}"

rule concat_gvcfs:
    input:
        get_gvcfs
    output:
        gvcf = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz",
        tbi = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz.tbi"
    conda:
        "../envs/qc.yml"
    shell:
        """
        bcftools concat -O z -o {output.gvcf} {input}
        tabix -p vcf {output.gvcf}
        """
    
