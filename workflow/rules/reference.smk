ruleorder: download_reference > index_reference
# localrules: copy_reference, download_reference

# This does not work with SLURM as of 4/3/24. See here for more info:https://github.com/snakemake/snakemake-executor-plugin-slurm/issues/60
# rule copy_reference:
#     """Copies user-specified reference genome path to results dir to maintain refGenome wildcard"""
#     input:
#         ref = get_ref
#     output:
#         ref = "results/{refGenome}/data/genome/{refGenome}.fna"
#     log:
#         "logs/{refGenome}/copy_ref/log.txt"
#     shell:
#         #probably don't need to unzip but might as well.
#         """
#         gunzip -c {input.ref} 2> {log} > {output.ref} || cp {input.ref} {output.ref} &> {log}
#         """

rule download_reference:
    input:
        ref = get_ref
    output:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna"
    params:
        dataset = "results/{refGenome}/data/genome/{refGenome}_dataset.zip",
        outdir = "results/{refGenome}/data/genome/{refGenome}"
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/download_ref/log.txt"
    benchmark:
        "benchmarks/{refGenome}/download_ref/benchmark.txt"
    shell:
        """
        if [ -z "{input.ref}" ]  # check if this is empty
        then
            mkdir -p {params.outdir}
            datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} \
            && 7z x {params.dataset} -aoa -o{params.outdir} \
            && cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}
        else
            gunzip -c {input.ref} 2> {log} > {output.ref} || cp {input.ref} {output.ref} &> {log}
        fi
        """

rule index_reference:
    input:
        ref = "results/{refGenome}/data/genome/{refGenome}.fna"
    output:
        indexes = expand("results/{{refGenome}}/data/genome/{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"]),
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
        dictf = "results/{refGenome}/data/genome/{refGenome}.dict"
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{refGenome}/index_ref/log.txt"
    benchmark:
        "benchmarks/{refGenome}/index_ref/benchmark.txt"
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} --output {output.fai} >> {log}
        samtools dict {input.ref} -o {output.dictf} >> {log} 2>&1
        """
