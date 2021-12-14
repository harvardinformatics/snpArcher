localrules: genome_prep

## RULES ##

rule genome_prep:
  input:
      ref = config["refGenomeDir"] + "{refGenome}.fna"
  output:
      twobit = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".2bit",
      chrom = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".sizes"
  conda:
      "../envs/callable.yml"
  resources:
      mem_mb = 
  shell:
      "faToTwoBit {input.ref} {output.twobit}\n"
      "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.sizes}"

rule bedgraphs:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix'],
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        bedgraph = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".bg"),
        sorted = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".sorted.bg")
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bedtools genomecov -ibam {input.bam} -bga -g {input.ref} > {output.bedgraph}\n"
        "sort k1,1 -k2,2n {output.bedgraph} > {output.sorted}"
        
rule merge:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna", 
        lambda wildcards:
        expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "{sample}.bam", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    output:
        merge = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".merge.bg")
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bedtools unionbedg -header -empty -g {input.ref} -names -i {input} > {output.merge}"
        
 rule bigBeds:
    input:
        merge = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".merge.bg",
        chrom = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".sizes"
    output:
        bigbed = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".merge.bigBed"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bedToBigBed {output.merge} {input.chrom} {output.bigbed}"

rule write_beds:
    input:
        bed = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".merge.bigBed"
    output:
        clean = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".clean.bed",
        high = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".high.bed",
        low = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".low.bed"
    conda:
        "../envs/callable.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "Rscript "
