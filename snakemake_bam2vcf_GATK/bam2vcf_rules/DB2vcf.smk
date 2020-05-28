rule DB2vcf:
    """
    This rule uses the genomic databases from the previous step (gvcf2DB) to create VCF files, one per list file. Thus, lists
    are still scattered.
    """
    input:
        DB = directory(dbDir + "DB_L{list}"),
        ref = config['ref']
    output: 
        vcf = vcfDir + "L{list}.vcf"
    resources: 
        mem_gb = int(CLUSTER["DB2vcf"]["mem"]/1000)
    run:
        command="""module load jdk/1.8.0_45-fasrc01
        /n/home11/bjarnold/gatk-4.1.0.0/gatk GenotypeGVCFs \
        --java-options \"-Xmx{resources.mem_gb}g -Xms{resources.mem_gb}g\" \
        -R {input.ref} \
        -V gendb://{input.DB} \
        -O {output.vcf} \
        --tmp-dir={vcfDir}tmp"""

        shell(command)
