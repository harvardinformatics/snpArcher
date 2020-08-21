rule processRef:
    """
    This rule generates a .fai file from the reference genome, which are required for GATK to run properly. GATK also needs a .dict file, but this was previously generated.
    """
    input:
        ref = config['ref'],
    output: 
        fai = config['ref'] + ".fai",
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = 5000
    shell:
        "samtools faidx {input.ref} --output {output.fai}\n"

rule bam2gvcf:
    """
    This rule scatters analyses over two dimensions: sample name and list file. For each BAM file, one per sample,
    a GVCF is created for all the scaffolds present in a given list file.
    """
    input:
        ref = config['ref'],
        fai = config['ref'] + ".fai",
        dict = refBaseName + ".dict",
        bam = bamDir + "{sample}_dedup.bam",
        l = listDir + "list{list}.list"
    output: 
        gvcf = gvcfDir + "{sample}_L{list}.raw.g.vcf.gz",
        gvcf_idx = gvcfDir + "{sample}_L{list}.raw.g.vcf.gz.tbi"
    resources: 
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['bam2gvcf']['mem'] - 3000)  # this is the maximum amount given to java
    params:
        minPrun = 1,
        minDang = 1
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk HaplotypeCaller "
        "--java-options \"-Xmx{resources.reduced}m\" "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-L {input.l} "
        "--emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang}"

rule gvcf2DB:
    """
    This rule gathers results for a given list file name, so the workflow is now scattered in only a single dimension. 
    Here we take many gvcfs for a particular list of scaffolds and combine them into a GenomicsDB.
    Samples are thus gathered by a shared list name, but lists are still scattered.
    """
    input:
        # NOTE: this waits for all gvcfs to be finished, whereas ou really only need to wait 
        # for all samples from a particular list to be finished
        gvcfs = expand(gvcfDir + "{sample}_L{list}.raw.g.vcf.gz", sample=SAMPLES, list=LISTS),
        gvcfs_idx = expand(gvcfDir + "{sample}_L{list}.raw.g.vcf.gz.tbi", sample=SAMPLES, list=LISTS),
        l = listDir + "list{list}.list",
        DBmapfile = dbDir + "DB_mapfile{list}"
    output: 
        DB = directory(dbDir + "DB_L{list}"),
        doneFile = temp(touch(dbDir + "DB_L{list}.done"))
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['gvcf2DB']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['gvcf2DB']['mem'] - 3000)  # this is the maximum amount given to java
    conda:
        "../envs/bam2vcf.yml"
    shell:
        # NOTE: reader-threads > 1 useless if you specify multiple intervals
        # a forum suggested TILEDB_DISABLE_FILE_LOCKING=1 to remedy sluggish performance
        "export TILEDB_DISABLE_FILE_LOCKING=1 \n"
        "gatk GenomicsDBImport "
        "--java-options \"-Xmx{resources.reduced}m -Xms{resources.reduced}m\" "
        "--genomicsdb-workspace-path {output.DB} "
        "-L {input.l} "
        "--tmp-dir {dbDir}tmp "
        "--sample-name-map {input.DBmapfile} \n"

rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        DB = dbDir + "DB_L{list}",
        ref = config['ref'],
        doneFile = dbDir + "DB_L{list}.done"
    output: 
        vcf = vcfDir + "L{list}.vcf"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['DB2vcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['DB2vcf']['mem'] - 3000)  # this is the maximum amount given to java
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk GenotypeGVCFs "
        "--java-options \"-Xmx{resources.reduced}m -Xms{resources.reduced}m\" "
        "-R {input.ref} "
        "-V gendb://{input.DB} "
        "-O {output.vcf} "
        "--tmp-dir {vcfDir}tmp"

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
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gatherVcfs']['mem']   # this is the overall memory requested
    shell:
        "INPUT=\"\" \n"
        "for ((i=0;i<={lastList};i++)) \n"
        "do \n"
        "INPUT=\"${{INPUT}} -I {vcfDir}L${{i}}.vcf\" \n"
        "done\n"

        "gatk GatherVcfs "
        "$INPUT "
        "-O {output.vcf}\n"

        "sleep 10\n" # the variant filtration step was failing in an unreproducible way, so added this in case
        # Hard filter Combined.vcf
        "gatk VariantFiltration "
        "-R {input.ref} "
        "-V {output.vcf} "
        "--output {output.vcfFiltered} "
        "--filter-expression \"QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || ExcessHet > 30.0\" "
        "--filter-name \"filteredOut\" "

rule vcftools:
    input:
        vcf = config["gatkDir"] + "Combined_hardFiltered.vcf",
        int = config["gatkDir"] + "intervals.bed"
    output: 
        missing = config["gatkDir"] + "missing_data_per_ind.txt",
        SNPsPerInt = config["gatkDir"] + "SNP_per_interval.txt"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        "vcftools --vcf {input.vcf} --remove-filtered-all --minDP 1 --stdout --missing-indv > {output.missing}\n"
        "bedtools intersect -a {input.int} -b {input.vcf} -c > {output.SNPsPerInt}"
