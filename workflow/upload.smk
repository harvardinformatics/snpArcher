import pandas as pd
from yaml import safe_load
localrules: all, setup_dirs, upload_vcfs, readme, upload_gvcfs, upload_qc, upload_bams
configfile: config["config_file"]


samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
with open(config["resource_config"], "r") as f:
    resources = safe_load(f)


def get_output(wc):
    if config["final_prefix"] == "":
        raise(WorkflowError("'final_prefix' is not set in config."))
    out = []
    genomes = samples['refGenome'].unique().tolist()
    sample_counts = samples.drop_duplicates(subset = ["BioSample", "refGenome"]).value_counts(subset=['refGenome'])  #get BioSample for each refGenome
    out.extend
    for ref in genomes:
        _samples = samples['BioSample'].unique().tolist()
        out.extend(expand("upload_done_files/{refGenome}/{prefix}/vcfs.txt", refGenome=ref, prefix=config['final_prefix']))
        out.extend(expand("upload_done_files/{refGenome}/{prefix}/qc_files.txt", refGenome=ref, prefix=config['final_prefix']))
        out.extend(expand("upload_done_files/{refGenome}/{prefix}/readme.txt", refGenome=ref, prefix=config['final_prefix']))
        out.extend(expand("upload_done_files/{refGenome}/{prefix}/gvcf_{sample}.txt", refGenome=ref, prefix=config['final_prefix'], sample=_samples))
        out.extend(expand("upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt", refGenome=ref, prefix=config['final_prefix']))
        # out.extend(expand("upload_done_files/{refGenome}/{prefix}/bam_{sample}.txt", refGenome=ref, prefix=config['final_prefix'], sample=_samples))
    
    return out

rule all:
    input: get_output

rule setup_dirs:
    output: touch("upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt")
    script: "scripts/upload.py"

rule upload_vcfs:
    input:
        "upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt",
        "upload_done_files/{refGenome}/{prefix}/qc_files.txt",
        "upload_done_files/{refGenome}/{prefix}/readme.txt",
        vcf = "results/{refGenome}/{prefix}_final.vcf.gz",
        tbi = "results/{refGenome}/{prefix}_final.vcf.gz.tbi",
        filtered = "results/{refGenome}/QC/{prefix}_filtered.vcf.gz",
        filtered_idx = "results/{refGenome}/QC/{prefix}_filtered.vcf.gz.csi",
    output:
        touch("upload_done_files/{refGenome}/{prefix}/vcfs.txt")
    script: "scripts/upload.py"

rule readme:
    input:
        "upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt",
    output:
        "results/{refGenome}/{prefix}_README.pdf",
        touch("upload_done_files/{refGenome}/{prefix}/readme.txt")
    script:
        "scripts/upload.py"

rule upload_gvcfs:
    input:
        # "upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt",
        # "upload_done_files/{refGenome}/{prefix}/qc_files.txt",
        # "upload_done_files/{refGenome}/{prefix}/readme.txt",
        # "upload_done_files/{refGenome}/{prefix}/vcfs.txt",
        "upload_done_files/{refGenome}/{prefix}/qc_files.txt",
        gvcf = "results/{refGenome}/gvcfs/{sample}.g.vcf.gz",
        gvcf_idx = "results/{refGenome}/gvcfs/{sample}.g.vcf.gz.tbi",
    output:
        touch("upload_done_files/{refGenome}/{prefix}/gvcf_{sample}.txt")
    script:
        "scripts/upload.py"

rule upload_bams:
    input:
        # "upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt",
        # "upload_done_files/{refGenome}/{prefix}/qc_files.txt",
        # "upload_done_files/{refGenome}/{prefix}/readme.txt",
        # "upload_done_files/{refGenome}/{prefix}/vcfs.txt",
        bam = "results/{refGenome}/bams/{sample}_final.bam",
        bai = "results/{refGenome}/bams/{sample}_final.bam.bai",
    output:
        touch("upload_done_files/{refGenome}/{prefix}/bam_{sample}.txt")
    script:
        "scripts/upload.py"

rule upload_qc:
    input:
        "upload_done_files/{refGenome}/{prefix}/readme.txt",
        "upload_done_files/{refGenome}/{prefix}/set_up_dirs.txt",
        eigenvec = "results/{refGenome}/QC/{prefix}.eigenvec",
        eigenval = "results/{refGenome}/QC/{prefix}.eigenval",
        depth = "results/{refGenome}/QC/{prefix}.idepth",
        dist = "results/{refGenome}/QC/{prefix}.dist",
        distid = "results/{refGenome}/QC/{prefix}.dist.id",
        king = "results/{refGenome}/QC/{prefix}.king",
        miss = "results/{refGenome}/QC/{prefix}.imiss",
        admix3 = "results/{refGenome}/QC/{prefix}.3.Q",
        admix2 = "results/{refGenome}/QC/{prefix}.2.Q",
        snpqc = "results/{refGenome}/QC/{prefix}_snpqc.txt",
        faiResult = "results/{refGenome}/QC/{prefix}_fai_tmp.txt",
        bed = "results/{refGenome}/QC/{prefix}.bed",
        bim = "results/{refGenome}/QC/{prefix}.bim",
        fam = "results/{refGenome}/QC/{prefix}.fam",
        sumstats = "results/{refGenome}/QC/{prefix}_bam_sumstats.txt",
        summ = "results/{refGenome}/QC/{prefix}.FILTER.summary",
        het = "results/{refGenome}/QC/{prefix}.het",
        fai = "results/{refGenome}/QC/{prefix}.fna.fai",
        qcpdf = "results/{refGenome}/QC/{prefix}_qc.html"
    output:
        touch("upload_done_files/{refGenome}/{prefix}/qc_files.txt"),
    script:
        "scripts/upload.py"





