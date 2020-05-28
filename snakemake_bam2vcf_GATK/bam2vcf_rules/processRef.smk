rule processRef:
    """
    This rule generates a .dict and .fai file from the reference genome, which are required for GATK to run properly.
    """
    input:
        ref = config['ref'],
    output: 
        fai = config['ref'] + ".fai",
        dict = refBasename + ".dict"
    run:
        shell("/n/home11/bjarnold/programs/samtools-1.10/samtools faidx {input.ref} --output {output.fai}") 
        shell("java -jar /n/home11/bjarnold/picard.jar CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output.dict}")
