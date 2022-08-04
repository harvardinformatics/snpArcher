rule check_fai:
    """
    checks fai file for numeric first column, then do not run plink and rest of workflow if they are all numeric
    """
    input:
        vcf = "results/{refGenome}/final.vcf.gz",
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
    output:
        faiResult = "results/{refGenome}/QC/{refGenome}_fai_tmp.txt"
    run:
        check_contig_names(input.fai, output.faiResult)

rule vcftools_individuals:
    input:
        vcf = "results/{refGenome}/final.vcf.gz"
    output:
        depth = "results/{refGenome}/QC/{refGenome}.idepth",
        miss = "results/{refGenome}/QC/{refGenome}.imiss",
        samps = "results/{refGenome}/QC/{refGenome}.samps.txt",
        summ = "results/{refGenome}/QC/{refGenome}.FILTER.summary",
        het = "results/{refGenome}/QC/{refGenome}.het"
    conda:
        "../envs/qc.yml"
    params:
        prefix = lambda wc, input: input.het[:-4],
        min_depth = config["min_depth"]
    shell:
        """
        vcftools --gzvcf {input.vcf} --FILTER-summary --out {params.prefix}
        vcftools --gzvcf {input.vcf} --out {params.prefix} --depth
        vcftools --gzvcf {input.vcf} --out {params.prefix} --het
        vcftools --gzvcf {input.vcf} --out {params.prefix} --missing-indv
        tail +1 {output.depth} | awk "$3>{params.min_depth}" {{print $1}}> {output.samps}
        """

rule subsample_snps:
    input:
        vcf = "results/{refGenome}/final.vcf.gz",
        samps = "results/{refGenome}/QC/{refGenome}.samps.txt",
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
        sumstats = "results/{refGenome}/summary_stats/bam_sumstats.txt"

    output:
        filtered = config['output'] + "{Organism}/{refGenome}/" + "{refGenome}.filtered.vcf.gz",
        filtered_idx = config['output'] + "{Organism}/{refGenome}/" + "{refGenome}.filtered.vcf.gz.csi",
        pruned = "results/{refGenome}/QC/{refGenome}.pruned.vcf.gz",
        snpqc = "results/{refGenome}/QC/{refGenome}_snpqc.txt",
        fai = "results/{refGenome}/QC/{refGenome}.fna.fai",
        sumstats = "results/{refGenome}/QC/{refGenome}_bam_sumstats.txt"
    conda:
        "../envs/qc.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['vcftools']['mem']
    shell:
        """
        ##first remove filtered sites and retain only biallelic SNPs
        ##Also remove sites with MAF < 0.01 and those with > 75% missing data
        bcftools view -S {input.samps} -t ^mtDNA -v snps -m2 -M2 -f .,PASS -e 'AF==1 | AF==0 | AF<0.01 | ALT="*" | F_MISSING > 0.75 | TYPE~"indel" | ref="N"' {input.vcf} -O z -o {output.filtered}
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

        bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\t%QUAL\t%INFO/ReadPosRankSum\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\n' {output.pruned} > {output.snpqc}
        
        ##copy the fai file into the QC folder for easy access
        cp {input.fai} {output.fai}
        cp {input.sumstats} {output.sumstats}
        """

rule plink:
    """
    Call plink PCA.
    """
    input:
        vcf = "results/{refGenome}/QC/{refGenome}.pruned.vcf.gz",
        faiResult = "results/{refGenome}/QC/{refGenome}_fai_tmp.txt"        
    params:
        prefix = os.path.join(workflow.default_remote_prefix, ("results/{refGenome}/QC/{refGenome}")),
        threads = config['plink']['threads']
    output: 
        bed = "results/{refGenome}/QC/{refGenome}.bed",
        bim = "results/{refGenome}/QC/{refGenome}.bim",
        fam = "results/{refGenome}/QC/{refGenome}.fam",
        eigenvec = "results/{refGenome}/QC/{refGenome}.eigenvec",
        eigenval = "results/{refGenome}/QC/{refGenome}.eigenval",
        dist = "results/{refGenome}/QC/{refGenome}.dist",
        distid = "results/{refGenome}/QC/{refGenome}.dist.id",
        king = "results/{refGenome}/QC/{refGenome}.king"

    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['vcftools']['mem']    # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    shell:
        #plink 2 for king relatedness matrix (robust to structure) and plink 1.9 for distance matrix
        """
        plink2 --vcf {input.vcf} --pca 10 --out {params.prefix} --allow-extra-chr --autosome-num 95 --make-bed --make-king square --const-fid --bad-freqs
        plink --vcf {input.vcf} --out {params.prefix} --allow-extra-chr --autosome-num 95 --distance square --const-fid
        """

rule admixture:
    """
    Call Admixture. First, make a bim file that has no charecters in the chromosomes
    """
    input:
        bed = "results/{refGenome}/QC/{refGenome}.bed",
        bim = "results/{refGenome}/QC/{refGenome}.bim",
        fam = "results/{refGenome}/QC/{refGenome}.fam",
    output:
        admix = "results/{refGenome}/QC/{refGenome}.3.Q",
        admix2 = "results/{refGenome}/QC/{refGenome}.2.Q"
    params:
        tmpbim = "results/{refGenome}/QC/{refGenome}",
        outdir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['qcDir']))
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['admixture']['mem'] 
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
        "results/{refGenome}/QC/{refGenome}.coords.txt"
    run:

        out_df = samples[["BioSample", "long", "lat"]]
        out_df.drop_duplicates("BioSample", inplace=True)
        out_df.dropna(subset=["long", "lat"], thresh=1, inplace=True)
        outpath = os.path.join(workflow.default_remote_prefix, config['output'], f"{wildcards.Organism}/{wildcards.ref}/", config['qcDir'], f"{wildcards.Organism}_{wildcards.refGenome}.coords.txt")
        out_df.to_csv(outpath, index=False, sep="\t")


rule qc_plots:
    """
    Call plotting script
    """
    input:
        eigenvec = "results/{refGenome}/QC/{refGenome}.eigenvec",
        eigenval = "results/{refGenome}/QC/{refGenome}.eigenval",
        depth = "results/{refGenome}/QC/{refGenome}.idepth",
        dist = "results/{refGenome}/QC/{refGenome}.dist",
        distid = "results/{refGenome}/QC/{refGenome}.dist.id",
        king = "results/{refGenome}/QC/{refGenome}.king",
        miss = "results/{refGenome}/QC/{refGenome}.imiss",
        admix3 = "results/{refGenome}/QC/{refGenome}.3.Q",
        admix2 = "results/{refGenome}/QC/{refGenome}.2.Q",
        snpqc = "results/{refGenome}/QC/{refGenome}_snpqc.txt",
        faiResult = "results/{refGenome}/QC/{refGenome}_fai_tmp.txt",
        bed = "results/{refGenome}/QC/{refGenome}.bed",
        bim = "results/{refGenome}/QC/{refGenome}.bim",
        fam = "results/{refGenome}/QC/{refGenome}.fam",
        sumstats = "results/{refGenome}/QC/{refGenome}_bam_sumstats.txt",
        summ = "results/{refGenome}/QC/{refGenome}.FILTER.summary",
        het = "results/{refGenome}/QC/{refGenome}.het",
        fai = "results/{refGenome}/QC/{refGenome}.fna.fai",
        coords = get_coords_if_available
    params:
        prefix = lambda wc, input: input.het[:-4],
        nClusters = config['nClusters'],
        GMKey = config['GoogleAPIKey']
    output: 
        qcpdf = "results/{refGenome}/QC/{refGenome}_qc.html"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['vcftools']['mem'] # this is the overall memory requested
    conda:
        "../envs/qc.yml"
    script:
        "../scripts/qc_dashboard_render.R"
