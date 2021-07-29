### INPUT FUNCTIONS ###
def get_read_group(wildcards):
    """Denote sample name and library_id in read group."""
    return r"-R '@RG\tID:{lib}\tSM:{sample}\tPL:ILLUMINA'".format(
        sample=wildcards.sample,
        lib=samples.loc[samples['BioSample'] == wildcards.sample]["LibraryName"].tolist()[0]
    )


def get_sumstats(wildcards):
    # Gets the correct sample given the organism and reference genome
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    fastpFiles = expand(config['output'] + config['sumstatDir'] + "{{Organism}}/{{refGenome}}/{sample}_fastp.out", sample=_samples)
    dedupFiles = expand(config['output'] + config['sumstatDir'] + "{{Organism}}/{{refGenome}}/{sample}_dedupMetrics.txt", sample=_samples)
    alnSumMetsFiles = expand(config['output'] + config['sumstatDir'] + "{{Organism}}/{{refGenome}}/{sample}_AlnSumMets.txt", sample=_samples)
    coverageFiles = expand(config['output'] + config['sumstatDir'] + "{{Organism}}/{{refGenome}}/{sample}_coverage.txt", sample=_samples)
    validateFiles = expand(config['output'] + config['sumstatDir'] + "{{Organism}}/{{refGenome}}/{sample}_validate.txt", sample=_samples)
    return {'fastpFiles': fastpFiles, 'dedupFiles': dedupFiles, 'alnSumMetsFiles': alnSumMetsFiles, 'coverageFiles': coverageFiles, 'validateFiles': validateFiles}