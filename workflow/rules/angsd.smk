rule bamlist:
    input: 
        bams = get_final_bams,
        vcf = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz"
    output:
        bamlist = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.bamlist"
    conda:
        "../envs/angsd.yml"
    shell:
        "ls {input.bams} > {output.bamlist}"

rule angsd:
    input:
        bamlist = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.bamlist",
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd.beagle.gz"
    params:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd"
    conda:
        "../envs/angsd.yml"
    threads: 31
    shell:
        """
        angsd -b {input.bamlist} -ref {input.ref}  \
                -out {params.angsd} \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
                -minMapQ 20 -minQ 20 \
                -GL 2 -P {threads} -doGlf 2 -doMajorMinor 4
        """

rule installpcangsd:
    input:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd.beagle.gz"
    output:
        confirm = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.confirm"
    conda:
        "../envs/angsd.yml"
    shell:
        """
        rm -r pcangsd
        git clone https://github.com/Rosemeis/pcangsd.git
        cd pcangsd
        pip3 install -e .
        cd ..
        echo "pcangsd install complete" > {output.confirm}
        """

rule pcangsd:
    input:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd.beagle.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.cov"
    params:
        cov = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.cov"
    conda:
        "../envs/angsd.yml"
    threads: 31
    shell:
        """
        pcangsd --beagle {input.angsd} \
        --out {output.cov} -t {threads}
        """