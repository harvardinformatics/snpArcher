localrules: genome_prep

## RULES ##

rule genome_prep:
  input:
      ref = config["refGenomeDir"] + "{refGenome}.fna"
  output:
      twobit = config['output'] + "{refGenome}/" + "{refGenome}" + ".2bit",
      chrom = config['output'] + "{refGenome}/" + "{refGenome}" + ".sizes"
  conda:
      "../envs/callable.yml"
  shell:
      "faToTwoBit {input.ref} {output.twobit}\n"
      "twoBitInfo {output.twobit} stdout | sort -k2rn > {output.chrom}"

rule bedgraphs:
    input:
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix']
    output:
        temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}.sorted.bg")
    conda:
        "../envs/callable.yml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bedtools genomecov -ibam {input.bam} -bga | sort -k1,1 -k2,2n - > {output}"

rule merge_bedgraph:
    input:
        unpack(get_input_for_coverage)
    output:
        merge = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{Organism}.merge.bg")
    log:
        "logs/{Organism}/covcalc/{refGenome}_{Organism}_merge.txt"
    conda:
        "../envs/callable.yml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bedtools unionbedg -header -empty -g {input.chrom} -i {input.bedgraphs} > {output.merge} 2> {log}"

rule gzip_bedgraph:
    input:
        get_bedgraph_to_convert
    output:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz",
        idx = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gzi"
    conda:
        "../envs/callable.yml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
    shell:
        "bgzip -i -I {output.idx} -c {input} > {output.bgz}"

#rule write_beds:
#    input:
#        bed = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".merge.bigBed"
#    output:
#        clean = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".clean.bed",
#        high = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".high.bed",
#        low = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + ".low.bed"
#    conda:
#        "../envs/callable.yml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * res_config['bedgraphs']['mem']
#    shell:
#        "Rscript "
