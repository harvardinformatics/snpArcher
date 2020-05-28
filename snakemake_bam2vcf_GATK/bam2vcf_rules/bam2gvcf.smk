rule bam2gvcf:
    """
    This rule scatters analyses over two dimensions: sample name and list file. For each BAM file, one per sample,
    a GVCF is created for all the scaffolds present in a given list file.
    """
    input:
        ref = config['ref'],
        fai = config['ref'] + ".fai",
        dict = refBasename + ".dict",
        bam = bamDir + "{sample}_dedupSort.bam",
        l = listDir + "list{list}.list"
    output: 
        gvcf = gvcfDir + "{sample}_L{list}.raw.g.vcf",
    resources: 
        cpus = CLUSTER["bam2gvcf"]["n"],
        mem_gb = int(CLUSTER["bam2gvcf"]["mem"]/1000)
    params:
        minPrun = 1,
        minDang = 1
    run:
        command = """module load jdk/1.8.0_45-fasrc01
        /n/home11/bjarnold/gatk-4.1.0.0/gatk HaplotypeCaller \
        --java-options \"-Xmx{resources.mem_gb}g -XX:ParallelGCThreads={resources.cpus}\" \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -L {input.l} \
        --emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang}"""

        shell(command)
