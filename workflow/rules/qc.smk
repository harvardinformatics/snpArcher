rule check_fai:
    """
    checks fai file for numeric first column, then do not run QC pipeline if they are all numeric
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai"
    output:
        faiResult = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_fai_tmp.txt"
    run:
        check_contig_names(input.fai, output.faiResult)

rule vcftools_individuals:
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    output:
        depth = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.idepth",
        miss = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.imiss"
    conda:
        "../envs/qc.yml"
    params:
        prefix= os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"))
    shell:
        """
        vcftools --gzvcf {input.vcf} --FILTER-summary --out {params.prefix}
        vcftools --gzvcf {input.vcf} --out {params.prefix} --depth
        vcftools --gzvcf {input.vcf} --out {params.prefix} --missing-indv
        """

rule subsample_snps:
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    output:
        filtered = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.filtered.vcf.gz",
        filtered_idx = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.filtered.vcf.gz.csi",
        pruned = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.pruned.vcf.gz"
    conda:
        "../envs/qc.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        """
        #first remove filtered sites and retain only biallelic variants
        bcftools view -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | ALT="*" | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
        bcftools index {output.filtered}
        #figure out how many SNPs are left, then identify how big of SNP window size to get down to between 100 and 150k snps        
        ALLSITES=`bcftools query -f '%CHROM\t%POS\n' {output.filtered} | wc -l`
        SITES=`echo $(( ${{ALLSITES}} / 100000 ))`

        #if the top VCF has < 150k SNPs, then just take all the SNPs
        if [[ $SITES -gt 1 ]]
        then
            bcftools +prune -w $SITES -n 1 -N rand -O z -o {output.pruned} {output.filtered}
        else
            bcftools view -O z -o {output.pruned} {output.filtered}
        fi
        """

rule snp_filters_qc:
    input:
        #I switched this input to be the pruned dataset, but this may be suboptimal for evaluating var quality distributions
        vcf = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.pruned.vcf.gz"
    output:
        snpqc = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_snpqc.txt"
    conda:
        "../envs/qc.yml"
    params:
        prefix=config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\t%QUAL\t%INFO/ReadPosRankSum\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\n' {input.vcf} > {output.snpqc}"


rule plink:
    """
    Call plink PCA.
    """
    input:
        vcf = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.pruned.vcf.gz",
        faiResult = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_fai_tmp.txt"        
    params:
        prefix = os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}")),
        threads = res_config['plink']['threads']
    output: 
        bed = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bed",
        bim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bim",
        fam = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.fam",
        eigenvec = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenvec",
        eigenval = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenval",
        dist = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.dist",
        distid = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.dist.id",
        king = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.king"

    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    shell:
        #plink 2 for king relatedness matrix (robust to structure) and plink 1.9 for distance matrix
        """
        plink2 --vcf {input.vcf} --pca 2 --out {params.prefix} --allow-extra-chr --autosome-num 95 --make-bed --make-king square --const-fid --bad-freqs
        plink --vcf {input.vcf} --out {params.prefix} --allow-extra-chr --autosome-num 95 --distance square --const-fid
        """

rule admixture:
    """
    Call Admixture. First, make a bim file that has no charecters in the chromosomes
    """
    input:
        bed = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bed",
        bim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bim",
        fam = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.fam",
    output:
        admix = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.3.Q",
        admix2 = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.2.Q"
    params:
        tmpbim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}",
        outdir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['qcDir']))
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['admixture']['mem'] 
    conda:
        "../envs/qc.yml"
    shell:
        """
        mv {input.bim} {input.bim}.orig
        paste <(cut -f 1 {input.bim}.orig | sed 's/[^0-9]//g') <(cut -f 2,3,4,5,6 {input.bim}.orig) >  {input.bim}

        admixture {input.bed} 2
        admixture {input.bed} 3

        mv "{wildcards.Organism}_{wildcards.refGenome}".2.* {params.outdir}
        mv "{wildcards.Organism}_{wildcards.refGenome}".3.* {params.outdir}
        """

rule generate_coords_file:
    output: 
        config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.coords.txt"
    run:
        write_coords_file(wildcards)


rule qc_plots:
    """
    Call plotting script
    """
    input:
        eigenvec = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenvec",
        eigenval = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.eigenval",
        depth = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.idepth",
        dist = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.dist",
        distid = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.dist.id",
        king = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.king",
        miss = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.imiss",
        admix3 = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.3.Q",
        admix2 = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.2.Q",
        snpqc = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_snpqc.txt",
        faiResult = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_fai_tmp.txt",
        bed = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bed",
        bim = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.bim",
        fam = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.fam",
        coords = get_coords_if_available
    params:
        prefix = os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}")),
        nClusters = config['nClusters']
    output: 
        qcpdf = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}_qc.html"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem'] # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/qc_dashboard_render.R"
