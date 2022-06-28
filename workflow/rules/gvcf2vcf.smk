localrules: mkDBmapfile

rule mkDBmapfile:
    """
    This rule makes DB map files for the GenomicsDBImport function. These files contain the names of all the samples for a particular
    list. This can be necessary because with many samples (thousands) the gatk command can get too long if you include all sample
    names on the command line, so long SLURM throws an error.
    """
    input:
        unpack(get_input_for_mapfile)
    output:
        dbMapFile = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_mapfile.txt"
    run:
        write_db_mapfile(wildcards)

rule gvcf2DB:
    """
    This rule gathers results for a given list file name, so the workflow is now scattered in only a single dimension.
    Here we take many gvcfs for a particular list of scaffolds and combine them into a GenomicsDB.
    Samples are thus gathered by a shared list name, but lists are still scattered.
    """
    input:
        l = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "db_intervals/{list}-scattered.interval_list",
        dbMapFile = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_mapfile.txt"
    output:
        db = temp(directory(config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}")),
        tar = temp(config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}.tar"),
    params:
        tmp_dir = config['tmp_dir']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gvcf2DB']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: int(attempt * res_config['gvcf2DB']['mem'] * 0.80) # this is the maximum amount given to java
    log:
        "logs/{Organism}/gvcf2DB/{refGenome}_{list}.txt"
    benchmark:
        "benchmarks/{Organism}/gvcf2DB/{refGenome}_{list}.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        # NOTE: reader-threads > 1 useless if you specify multiple intervals
        # a forum suggested TILEDB_DISABLE_FILE_LOCKING=1 to remedy sluggish performance
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '-Xmx{resources.reduced}m -Xms{resources.reduced}m' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            -L {input.l} \
            --tmp-dir {params.tmp_dir} \
            --sample-name-map {input.dbMapFile} &> {log}
        
        tar --overwrite -cf {output.tar} {output.db}
        """

rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        db = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}.tar",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
    output:
        vcf = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf"),
        vcfidx = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf.idx")
    params:
        tmp_dir = config['tmp_dir'],
        het = config['het_prior'],
        db = config['output'] + "{Organism}/{refGenome}/" + config['dbDir'] + "DB_L{list}"
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
        """
        tar --overwrite -xf {input.db}

        gatk GenotypeGVCFs \
            --java-options '-Xmx{resources.reduced}m -Xms{resources.reduced}m' \
            -R {input.ref} \
            --heterozygosity {params.het} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {params.tmp_dir} &> {log}
        """

rule filterVcfs:
    """
    This rule filters all of the VCFs
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "L{list}.vcf",
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        vcf = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "filtered_L{list}.vcf.gz"),
        vcfidx = temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "filtered_L{list}.vcf.gz.tbi")
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
        "--create-output-variant-index  "
        "--invalidate-previous-filters true &> {log}"
rule make_list_of_vcfs:
    """Creates a file with list of vcfs for sort_gatherVcfs"""
    input:
        unpack(get_gather_vcfs)
    output:
        temp(config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "list_of_vcfs.txt")
    run:
        with open(output[0], 'w') as f:
            for line in input['gvcfs']:
                print(line, file=f)

rule sort_gatherVcfs:
    input:
        unpack(get_gather_vcfs),
        fof = config['output'] + "{Organism}/{refGenome}/" + config["vcfDir_gatk"] + "list_of_vcfs.txt",
    output:
        vcfFinal = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/{Organism}/sortVcf/{refGenome}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gatherVcfs']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['gatherVcfs']['mem'] - 2000)  # this is the maximum amount given to java
    benchmark:
        "benchmarks/{Organism}/sortVcf/{refGenome}.txt"
    shell:
        """
        bcftools concat -f {input.fof} -D -a -Ou | bcftools sort -Oz -o {output} -
        tabix -p vcf {output}
        """
