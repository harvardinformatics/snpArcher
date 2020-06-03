rule processRef:
    """
    This rule generates a .dict and .fai file from the reference genome, which are required for GATK to run properly.
    """
    input:
        ref = config['ref'],
    output: 
        fai = config['ref'] + ".fai",
        dict = refBaseName + ".dict"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "samtools faidx {input.ref} --output {output.fai}\n"
        "picard CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.dict}"

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
        gvcf = gvcfDir + "{sample}_L{list}.raw.g.vcf",
    resources: 
        cpus = CLUSTER["bam2gvcf"]["n"],
        mem_gb = int(CLUSTER["bam2gvcf"]["mem"]/1000)
    params:
        minPrun = 1,
        minDang = 1
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk HaplotypeCaller "
        "--java-options \"-Xmx{resources.mem_gb}g -XX:ParallelGCThreads={resources.cpus}\" "
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
        gvcfs = expand(gvcfDir + "{sample}_L{list}.raw.g.vcf", sample=SAMPLES, list=LISTS),
        l = listDir + "list{list}.list",
        DBmapfile = dbDir + "DB_mapfile{list}"
    output: 
        DB = directory(dbDir + "DB_L{list}")
    resources: 
        mem_gb = int(CLUSTER["gvcf2DB"]["mem"]/1000)
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk GenomicsDBImport "
        "--java-options \"-Xmx{resources.mem_gb}g -Xms{resources.mem_gb}g\" "
        "--genomicsdb-workspace-path {output.DB} "
        "-L {input.l} "
        "--tmp-dir={dbDir}tmp "
        "--sample-name-map {input.DBmapfile}"

rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        #DB = directory(dbDir + "DB_L{list}"),
        DB = dbDir + "DB_L{list}",
        ref = config['ref']
    output: 
        vcf = vcfDir + "L{list}.vcf"
    resources: 
        mem_gb = int(CLUSTER["DB2vcf"]["mem"]/1000)
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk GenotypeGVCFs "
        "--java-options \"-Xmx{resources.mem_gb}g -Xms{resources.mem_gb}g\" "
        "-R {input.ref} "
        "-V gendb://{input.DB} "
        "-O {output.vcf} "
        "--tmp-dir={vcfDir}tmp"

rule gatherVcfs:
    """
    This rule gathers all of the VCFs, one per list, into one final VCF
    """
    input:
        vcfs = expand(vcfDir + "L{list}.vcf", list=LISTS),
        ref = config['ref']
    output: 
        vcf = "Combined.vcf",
        vcfFiltered = "Combined_hardFiltered.vcf"
    resources: 
        mem_gb = int(CLUSTER["gatherVcfs"]["mem"]/1000)
    conda:
        "../envs/bam2vcf.yml"
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
