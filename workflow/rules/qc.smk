rule snp_filters_qc:
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    output:
        snpqc = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_snpqc.txt"
    conda:
        "../envs/qc.yml"
    params:
        prefix=config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\t%QUAL\t%INFO/ReadPosRankSum\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\n' {input.vcf} > {output.snpqc}"

rule vcftools_individuals:
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz",
        #vcf = "data/vcf_example/multisample_vcf.vcf.gz"
    output:
        depth = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.idepth",
        miss = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.imiss"
    conda:
        "../envs/qc.yml"
    params:
        prefix=config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"
    shell:
        """
        vcftools --gzvcf {input.vcf} --out {params.prefix} --depth
        vcftools --gzvcf {input.vcf} --out {params.prefix} --missing-indv
        """

rule subsample_snps:
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz",
    output:
        pruned = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.pruned.vcf.gz",
    conda:
        "../envs/qc.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        "bcftools +prune -n 1 -w 10000bp --nsites-per-win-mode rand --random-seed 42 -O v -o {output.pruned} {input.vcf}"

rule plink:
    """
    Call plink PCA.
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
        #vcf = "data/vcf_example/multisample_vcf.vcf.gz"
    params:
        prefix=config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"
    output: 
        bed = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bed",
        bim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bim",
        eigenvec = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenvec"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    shell:
        #plink 2 for king relatedness matrix (robust to structure) and plink 1.9 for distance matrix
        """
        plink2 --vcf {input.vcf} --pca 10 --out {params.prefix} --allow-extra-chr --autosome-num 30 --make-bed --make-king square --const-fid
        plink --vcf {input.vcf} --out {params.prefix} --allow-extra-chr --autosome-num 30 --distance square --const-fid
        """

rule admixture:
    """
    Call Admixture. First, make a bim file that has no charecters in the chromosomes
    """
    input:
        bed = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bed",
        bim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bim"
    output:
        admix = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.3.Q"
    params:
        tmpbim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}",
        outdir = config['output'] + "{Organism}/{refGenome}/" + config['qcDir']

    conda:
        "../envs/qc.yml"
    shell:
        """
        cut -f 1 {input.bim} > {params.tmpbim}_contig_names.tmp #get names of contigs
        sed 's/[^0-9]//g' {params.tmpbim}_contig_names.tmp > {params.tmpbim}_new_names.tmp #remove any letters from their names
        cp {input.bim} {params.tmpbim}.bim.bak #make a backup of the original bim (probably not neccessary)
        cut -f 2,3,4,5,6 {input.bim} > {params.tmpbim}_data.tmp # grab other columns from the bim file
        paste {params.tmpbim}_new_names.tmp {params.tmpbim}_data.tmp > {input.bim} # add updated names (no charecters) to the new bim, give it same name as original bim

        rm {params.tmpbim}_contig_names.tmp
        rm {params.tmpbim}_new_names.tmp
        rm {params.tmpbim}_data.tmp

        admixture {input.bed} 2
        admixture {input.bed} 3

        mv "{wildcards.Organism}_{wildcards.refGenome}".2.* {params.outdir}
        mv "{wildcards.Organism}_{wildcards.refGenome}".3.* {params.outdir}
        """

rule qc_plots:
    """
    Call plotting script
    """
    input:
        eigenvec = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenvec",
        depth = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.idepth",
        admix = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.3.Q",
        snpqc = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_snpqc.txt",
    params:
        prefix = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"
    output: 
        qcpdf = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_qc.html"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem'] # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/qc_dashboard_render.R"
