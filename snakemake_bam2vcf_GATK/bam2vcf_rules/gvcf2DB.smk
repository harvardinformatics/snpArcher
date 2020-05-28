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
    run:
        command="""module load jdk/1.8.0_45-fasrc01
        /n/home11/bjarnold/gatk-4.1.0.0/gatk GenomicsDBImport \
        --java-options \"-Xmx{resources.mem_gb}g -Xms{resources.mem_gb}g\" \
        --genomicsdb-workspace-path {output.DB} \
        -L {input.l} \
        --tmp-dir={dbDir}tmp \
        --sample-name-map {input.DBmapfile}"""

        shell(command)
