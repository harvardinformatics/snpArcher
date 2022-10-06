ruleorder: download_reference > index_reference

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
            cp {input.ref} {output.ref}
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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['index_ref']['mem']
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
