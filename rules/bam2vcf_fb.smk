
rule intervals:
    """
    Because snakemake input MUST be a file, this hack is necessary to iterate through intervals in the next step,
    inputing one interval in the command at a time. This information is stored in intervals_fb
    """
    output:
        interval = temp(touch(intervalDir + "{i}"))

rule bam2vcf:
    """
    Call freebayes.
    """
    input:
        ref = config['ref'],
        bams = expand(bamDir + "{sample}" + bam_suffix, sample=SAMPLES) # input all bams at once
    output: 
        vcf = temp(vcfDir_fb + "{i}.vcf") # output one vcf per interval
    conda:
        "../envs/bam2vcf.yml"
    resources:
        # increase memory every attempt
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2vcf']['mem'] 
    shell:
        "freebayes -f {input.ref} "
        "-r {wildcards.i} "
        "-g {maxDP_fb} "
        "{input.bams} "
        "> {output.vcf}\n"


rule gatherVcfs:
    """
    This rule gathers all of the VCFs, one per list, into one final VCF
    """
    input:
        vcfs = expand(vcfDir + "L{list}.vcf", list=LISTS),
        ref = config['ref']
    output: 
        vcf =  temp(config["gatkDir"] + "Combined.vcf"),
        vcfidx =  temp(config["gatkDir"] + "Combined.vcf.idx"),
        vcfFiltered =  config["gatkDir"] + "Combined_hardFiltered.vcf"
    params:
        gatherVcfsInput = helperFun.getVcfs_gatk(LISTS, vcfDir)
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gatherVcfs']['mem']   # this is the overall memory requested
    shell:
        "gatk GatherVcfs "
        "{params.gatherVcfsInput} "
        "-O {output.vcf}\n"

        "gatk VariantFiltration "
        "-R {input.ref} "
        "-V {output.vcf} "
        "--output {output.vcfFiltered} "
        "--filter-name \"RPRS_filter\" "
        "--filter-expression \"(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)\" "
        "--filter-name \"FS_SOR_filter\" "
        "--filter-expression \"(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))\" "
        "--filter-name \"MQ_filter\" "
        "--filter-expression \"vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))\" "
        "--invalidate-previous-filters"

rule compress:
    input:
        vcf = config["gatkDir"] + "Combined_hardFiltered.vcf"
    output:
        vcf = config["gatkDir"] + "Combined_hardFiltered.vcf.gz",
        idx = config["gatkDir"] + "Combined_hardFiltered.vcf.gz.tbi"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        "bgzip {input.vcf} && tabix -p vcf {output.vcf}"

rule vcftools:
    input:
        vcf = config["gatkDir"] + "Combined_hardFiltered.vcf.gz",
        int = intDir + "intervals_fb.bed"
    output: 
        missing = gatkDir + "missing_data_per_ind.txt",
        SNPsPerInt = gatkDir + "SNP_per_interval.txt"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        "vcftools --gzvcf {input.vcf} --remove-filtered-all --minDP 1 --stdout --missing-indv > {output.missing}\n"
        "bedtools intersect -a {input.int} -b {input.vcf} -c > {output.SNPsPerInt}"
