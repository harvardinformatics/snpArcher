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
    run:
        # combine VCFs	
        vcfList = "" 
        for i in input.vcfs:
                vcfList = vcfList + "-I {} ".format(i)
        command1="""module load jdk/1.8.0_45-fasrc01
        /n/home11/bjarnold/gatk-4.1.0.0/gatk GatherVcfs \
        {vcfList} \
        -O {output.vcf}"""

        shell(command1)
        shell("sleep 10") # the variant filtration step was failing in an unreproducible way, so added this in case
        # Hard filter Combined.vcf
        command2="""module load jdk/1.8.0_45-fasrc01
        /n/home11/bjarnold/gatk-4.1.0.0/gatk VariantFiltration \
        -R {input.ref} \
        -V {output.vcf} \
        --output {output.vcfFiltered} \
        --filter-expression \"QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || ExcessHet > 30.0\" \
        --filter-name \"filteredOut\"
        """

        shell(command2)
