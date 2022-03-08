localrules: mkDBmapfile

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
        int = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed"
    output:
        gvcf = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz",
        gvcf_idx = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz.tbi",
        doneFile = touch(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.done")
    resources:
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['bam2gvcf']['mem'] - 3000)  # this is the maximum amount given to java
    log:
        "logs/{Organism}/bam2gvcf/{refGenome}_{sample}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/bam2gvcf/{refGenome}_{sample}_{list}.txt"
    params:
        minPrun = config['minP'],
        minDang = config['minD'],
        l = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "list{list}.list"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk HaplotypeCaller "
        "--java-options \"-Xmx{resources.reduced}m\" "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-L {params.l} "
        "--emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang} &> {log}"

rule mkDBmapfile:
    """
    This rule makes DB map files for the GenomicsDBImport function. These files contain the names of all the samples for a particular
    list. This can be necessary because with many samples (thousands) the gatk command can get too long if you include all sample
    names on the command line, so long SLURM throws an error.
    """
    input:
        # NOTE: the double curly brackets around 'list' prevent the expand function from operating on that variable
        # thus, we expand by sample but not by list, such that we gather by sample for each list value
        unpack(get_input_for_mapfile)
    output:
        dbMapFile = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_mapfile_L{list}"
    run:
        write_db_mapfile(wildcards)

rule gvcf2DB:
    """
    This rule gathers results for a given list file name, so the workflow is now scattered in only a single dimension.
    Here we take many gvcfs for a particular list of scaffolds and combine them into a GenomicsDB.
    Samples are thus gathered by a shared list name, but lists are still scattered.
    """
    input:
        l = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "list{list}.list",
        dbMapFile = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_mapfile_L{list}"
    output:
        DB = directory(config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}"),
        done = touch(config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}.done")
    params:
        tmp_dir = config['tmp_dir']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gvcf2DB']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['gvcf2DB']['mem'] - 3000)  # this is the maximum amount given to java
    log:
        "logs/{Organism}/gvcf2DB/{refGenome}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/gvcf2DB/{refGenome}_{list}.txt"
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
        "--tmp-dir {params.tmp_dir} "
        "--sample-name-map {input.dbMapFile} &> {log}"

rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        DB = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        doneFile = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}.done"
    output:
        vcf = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf"),
        vcfidx = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf.idx")
    params:
        tmp_dir = config['tmp_dir']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['DB2vcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['DB2vcf']['mem'] - 3000)  # this is the maximum amount given to java
    log:
        "logs/{Organism}/db2vcf/{refGenome}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/db2vcf/{refGenome}_{list}.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        "gatk GenotypeGVCFs "
        "--java-options \"-Xmx{resources.reduced}m -Xms{resources.reduced}m\" "
        "-R {input.ref} "
        "-V gendb://{input.DB} "
        "-O {output.vcf} "
        "--tmp-dir {params.tmp_dir} &> {log}"

rule filterVcfs:
    """
    This rule filters all of the VCFs
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf",
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        vcf = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "filtered_L{list}.vcf"),
        vcfidx = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "filtered_L{list}.vcf.idx")
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['filterVcfs']['mem']   # this is the overall memory requested
    log:
        "logs/{Organism}/filterVcfs/{refGenome}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/filterVcfs/{refGenome}_{list}.txt"
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

rule sort_gatherVcfs:
    input: 
        get_gather_vcfs
    output: 
        vcfFinal = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    params:
        gather_vcfs_CLI
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gatherVcfs']['mem']
    log:
        "logs/{Organism}/sortVcf/{refGenome}.txt"
    benchmark:
        "benchmarks/{Organism}/sortVcf/{refGenome}.txt"
    shell:
        "gatk SortVcf "
        "{params} "
        "-O {output.vcfFinal} &> {log}"
