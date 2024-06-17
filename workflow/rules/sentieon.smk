rule sentieon_map:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        r1 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_1.fastq.gz",
        r2 = "results/{refGenome}/filtered_fastqs/{sample}/{run}_2.fastq.gz",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"])
    output: 
        bam = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam"),
        bai = temp("results/{refGenome}/bams/preMerge/{sample}/{run}.bam.bai"),
    params:
        rg = get_read_group,
        lic = config['sentieon_lic']
    conda:
        "../envs/sentieon.yml"
    log:
        "logs/{refGenome}/sentieon_map/{sample}/{run}.txt"
    benchmark:
        "benchmarks/{refGenome}/sentieon_map/{sample}/{run}.txt"
    shell:
        """
        export MALLOC_CONF=lg_dirty_mult:-1
        export SENTIEON_LICENSE={params.lic}
        sentieon bwa mem -M -R {params.rg} -t {threads} -K 10000000 {input.ref} {input.r1} {input.r2} | sentieon util sort --bam_compression 1 -r {input.ref} -o {output.bam} -t {threads} --sam2bam -i -
        samtools index {output.bam} {output.bai}
        """
rule merge_bams:
    input:
        merge_bams_input
    output:
        bam = temp("results/{refGenome}/bams/postMerge/{sample}.bam"),
        bai = temp("results/{refGenome}/bams/postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/merge_bams/{sample}.txt"
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam}"

rule sentieon_dedup:
    input:
        unpack(dedup_input),
    output:
        dedupBam = "results/{refGenome}/bams/{sample}_final.bam",
        dedupBai = "results/{refGenome}/bams/{sample}_final.bam.bai",
        score = temp("results/{refGenome}/summary_stats/{sample}/sentieon_dedup_score.txt"),
        metrics = temp("results/{refGenome}/summary_stats/{sample}/sentieon_dedup_metrics.txt")
    params:
        lic = config['sentieon_lic']
    conda:
        "../envs/sentieon.yml"
    log:
        "logs/{refGenome}/sentieon_dedup/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/sentieon_dedup/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.score}
        sentieon driver -t {threads} -i {input.bam} --algo Dedup --score_info {output.score} --metrics {output.metrics} --bam_compression 1 {output.dedupBam}
        """

rule sentieon_haplotyper:
    input:
        unpack(get_bams),
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        dictf = "results/{refGenome}/data/genome/{refGenome}.dict",
    params:
        lic = config['sentieon_lic'],
        ploidy = config['ploidy']
    output:
        gvcf = "results/{refGenome}/gvcfs/{sample}.g.vcf.gz",
        gvcf_idx = "results/{refGenome}/gvcfs/{sample}.g.vcf.gz.tbi",
    conda:
        "../envs/sentieon.yml"
    log:
        "logs/{refGenome}/sentieon_haplotyper/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/sentieon_haplotyper/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 {output.gvcf} --ploidy {params.ploidy} 2> {log}
        """

rule sentieon_combine_gvcf:
    input:
        unpack(sentieon_combine_gvcf_input),
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        dictf = "results/{refGenome}/data/genome/{refGenome}.dict"
    output:
        vcf = temp("results/{refGenome}/vcfs/raw.vcf.gz"),
        tbi = temp("results/{refGenome}/vcfs/raw.vcf.gz.tbi")
    params:
        glist = lambda wc, input: " ".join(["-v " + gvcf for gvcf in input['gvcfs']]),
        lic = config['sentieon_lic']
    conda:
        "../envs/sentieon.yml"
    log:
        "logs/{refGenome}/sentieon_combine_gvcf/log.txt"
    benchmark:
        "benchmarks/{refGenome}/sentieon_combine_gvcf/benchmark.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} -t {threads} --algo GVCFtyper --emit_mode VARIANT {output.vcf} {params.glist} 2> {log}
        """

rule filter_vcf:
    """
    This rule applies filters to the raw vcf using GNU Parallel.
    """
    input:
        vcf = "results/{refGenome}/vcfs/raw.vcf.gz",
        tbi = "results/{refGenome}/vcfs/raw.vcf.gz.tbi",
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        dictf = "results/{refGenome}/data/genome/{refGenome}.dict"
    output:
        vcf = "results/{refGenome}/{prefix}_raw.vcf.gz",
        tbi = "results/{refGenome}/{prefix}_raw.vcf.gz.tbi"
    conda:
        "../envs/bam2vcf.yml"
    log:
        "logs/{refGenome}/sentieon_filter_vcfs/{prefix}_log.txt"
    shadow: "minimal"
    benchmark:
        "benchmarks/{refGenome}/sentieon_filter_vcfs/{prefix}_benchmark.txt"
    shell:
        """
        # get the contig names from the .fai index
        contigs=$(cut -f1 {input.indexes[5]})
        
        # create a function that will be passed to gnu parallel
        filter_contig() {{
            contig=$1
            echo $contig

            gatk --java-options "-Xmx4g" VariantFiltration \
                -R {input.ref} \
                -L ${{contig}} \
                -V {input.vcf} \
                --output {wildcards.refGenome}_{wildcards.prefix}_filter_${{contig}}.vcf.gz \
                --filter-name "RPRS_filter" \
                --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" \
                --filter-name "FS_SOR_filter" \
                --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \
                --filter-name "MQ_filter" \
                --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
                --filter-name "QUAL_filter" \
                --filter-expression "QUAL < 30.0" \
                --invalidate-previous-filters true
        }}
        
        export -f filter_contig
        
        # pass each contig to gnu parallel
        parallel -j {threads} filter_contig ::: ${{contigs}} 2> {log}
        
        bcftools concat {wildcards.refGenome}_{wildcards.prefix}_filter_*.vcf.gz --threads {threads} -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """
