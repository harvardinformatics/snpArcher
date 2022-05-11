rule angsd:
    input:
        pruned = config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.pruned.vcf.gz",
        bams = get_final_bams,
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd.beagle.gz",
        bamlist = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.bamlist"
    params:
        angsd = os.path.join(workflow.default_remote_prefix, (config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd"))
    conda:
        "../envs/angsd.yml"
    threads: 31
    resources:
        disk_mb = 2048000,
        mem_mb = 512000
    shell:
        """
        samtools faidx {input.ref}
        ls -1 {input.bams} > {output.bamlist}

        angsd -b {output.bamlist} -ref {input.ref}  \
                -out {params.angsd} \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
                -minMapQ 20 -minQ 20 \
                -GL 2 -P {threads} -doGlf 2 -doMajorMinor 4
        """

rule pcangsd:
    input:
        angsd = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}_angsd.beagle.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.cov",
        confirm = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.confirm",
        args = config['output'] + "{Organism}/{refGenome}/" + config['angsdDir'] + "{Organism}_{refGenome}.cov.args"
    conda:
        "../envs/angsd.yml"
    threads: 31
    shell:
        """
        if [ -d pcangsd ]
        then
            echo "pcangsd install complete" > {output.confirm}
        else
            git clone https://github.com/Rosemeis/pcangsd.git
            cd pcangsd
            pip3 install -e .
            cd ..
            echo "pcangsd install complete" > {output.confirm}
            cd ..
        fi

        pcangsd --beagle {input.angsd} \
        --out {output.cov} -t {threads}
        """