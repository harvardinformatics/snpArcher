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

rule map:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        r1 = config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_1.fastq.gz",
        r2 = config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}_2.fastq.gz",
        indices = ancient(expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["fai", "sa", "pac", "bwt", "ann", "amb"])),
        lic = ancient(config['sentieon_lic'])
    output: 
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}/{run}.bam"),
        bai = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}/{run}.bam.bai"),
    params:
        rg = get_read_group_sentieon,
    conda:
        "../envs/sentieon.yml"
    threads: res_config['bwa_map']['threads']
    log:
        "logs/{Organism}/{refGenome}/bwa/{sample}_{run}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bwa_map']['mem'],
        machine_type = "n2d-standard-32"
    benchmark:
        "benchmarks/{Organism}/{refGenome}/bwa/{sample}/{run}.txt"
    shell:
        """
        export MALLOC_CONF=lg_dirty_mult:-1
        export SENTIEON_LICENSE={input.lic}
        sentieon bwa mem -M -R {params.rg} -t {threads} -K 10000000 {input.ref} {input.r1} {input.r2} 2>{log} | sentieon util sort --bam_compression 1 -r {input.ref} -o {output.bam} -t {threads} --sam2bam -i - 2>>{log}
        samtools index {output.bam} {output.bai} 2>>{log}
        """

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
        bam = get_bams_for_dedup,
        bai = get_bai_for_dedup,
        lic = ancient(config['sentieon_lic'])
    output:
        dedupBam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam"),
        dedupBai = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam.bai"),
    conda:
        "../envs/sentieon.yml"
    threads: 
        res_config['dedup']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['dedup']['mem'],
        machine_type = "n2d-standard-32"
    benchmark:
        "benchmarks/{Organism}/{refGenome}/dedup/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={input.lic}
        sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info score.txt
        sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt  --bam_compression 1 {output.dedupBam}
        """

rule gvcf:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        indices = ancient(expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["fai", "sa", "pac", "bwt", "ann", "amb"])),
        dictf = ancient(config["refGenomeDir"] + "{refGenome}" + ".dict"),
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam",
        bai = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + "_final.bam.bai",
        lic = ancient(config['sentieon_lic'])
    output:
        gvcf = temp(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz"),
        gvcf_idx = temp(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz.tbi"),
    threads: 31
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],
        machine_type = "n2d-standard-32"
    conda:
        "../envs/sentieon.yml"
    benchmark:
        "benchmarks/{Organism}/{refGenome}/gvcf/{sample}.txt"
    log:
        "logs/{Organism}/{refGenome}/bam2gvcf/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={input.lic}
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 {output.gvcf} 2> {log}
        """

rule combine_gvcf:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        indices = ancient(expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["fai", "sa", "pac", "bwt", "ann", "amb"])),
        dictf = ancient(config["refGenomeDir"] + "{refGenome}" + ".dict"),
        lic = ancient(config['sentieon_lic']),
        gvcfs = get_gvcfs,
        tbis = get_tbis
    output:
        vcf = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_prefilter.vcf.gz"),
        tbi = temp(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_prefilter.vcf.gz.tbi"),
    params:
        gvcf = get_gvcf_cmd
    threads: 31
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],
        machine_type = "n2d-standard-32"
    conda:
        "../envs/sentieon.yml"
    log: "logs/{Organism}/{refGenome}/combine_gvcf.txt"
    benchmark:
        "benchmarks/{Organism}/{refGenome}/combine_gvcf/{Organism}_{refGenome}.final.txt"
    shell:
        """
        export SENTIEON_LICENSE={input.lic}
        sentieon driver -r {input.ref} -t {threads} --algo GVCFtyper --emit_mode VARIANT {output.vcf} {params.gvcf} 2> {log}
        """

rule filterVcfs:
    """
    This rule filters all of the VCFs
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_prefilter.vcf.gz",
        tbi = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_prefilter.vcf.gz.tbi",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        indices = ancient(expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["fai", "sa", "pac", "bwt", "ann", "amb"])),
        dictf = ancient(config["refGenomeDir"] + "{refGenome}" + ".dict"),
    output:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['filterVcfs']['mem'],   # this is the overall memory requested
        machine_type = "n2d-standard-32"
    log:
        "logs/{Organism}/{refGenome}/filterVcfs/log.txt"
    benchmark:
        "benchmarks/{Organism}/{refGenome}/filterVcfs/bench.txt"
    shell:
        "gatk VariantFiltration "
        "-R {input.ref} "
        "-V {input.vcf} "
        "--output {output.vcf} "
        "--filter-name \"RPRS_filter\" "
        "--filter-expression \"(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)\" "
        "--filter-name \"FS_SOR_filter\" "
        "--filter-expression \"(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))\" "
        "--filter-name \"MQ_filter\" "
        "--filter-expression \"vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))\" "
        "--filter-name \"QUAL_filter\" "
        "--filter-expression \"QUAL < 30.0\" "
        "--invalidate-previous-filters true &> {log}"