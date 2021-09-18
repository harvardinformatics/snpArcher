import glob
import re
import os
from collections import defaultdict, deque

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
    fastpFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_fastp.out", sample=_samples)
    dedupFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_dedupMetrics.txt", sample=_samples)
    alnSumMetsFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_AlnSumMets.txt", sample=_samples)
    coverageFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_coverage.txt", sample=_samples)
    validateFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_validate.txt", sample=_samples)
    return {'fastpFiles': fastpFiles, 'dedupFiles': dedupFiles, 'alnSumMetsFiles': alnSumMetsFiles, 'coverageFiles': coverageFiles, 'validateFiles': validateFiles}

def get_gather_vcfs(wildcards):
    """
    Gets filtered vcfs for gathering step. This function gets the interval list indicies from the corresponding
    genome, then produces the file names for the filtered vcf with list index."""
    checkpoint_output = checkpoints.create_intervals.get(**wildcards).output[0]
    list_dir_search = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['intDir'], "*.list")
    list_files = glob.glob(list_dir_search)
    out = []
    for f in list_files:
        f = os.path.basename(f)
        index = re.search("\d+", f).group() # Grab digits from list file name and put in out list
        vcf = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['vcfDir_gatk'], f"filtered_L{index}.vcf")
        out.append(vcf)

    return out

def gather_vcfs_CLI(wildcards):
    """
    Gatk enforces that you have a -I before each input vcf, so this function makes that string
    """
    vcfs = get_gather_vcfs(wildcards)
    out = " ".join(["-I " + vcf for vcf in vcfs])
    return out

def write_db_mapfile(wildcards):
    dbMapFile = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['dbDir'], f"DB_mapfile_L{wildcards.list}")
    sample_names = set(samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist())
    with open(dbMapFile, 'w') as f:
        for sample in sample_names:
            gvcf_path = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['gvcfDir'], sample, f"L{wildcards.list}.raw.g.vcf.gz") 
            print(sample, gvcf_path, sep="\t", file=f)

def get_input_for_mapfile(wildcards):
    sample_names = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    gvcfs = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz", sample=sample_names, **wildcards)
    gvcfs_idx = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz.tbi", sample=sample_names, **wildcards)
    doneFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.done", sample=sample_names, **wildcards)
    print({'gvcfs': gvcfs, 'gvcfs_idx': gvcfs_idx, 'doneFiles': doneFiles})
    return {'gvcfs': gvcfs, 'gvcfs_idx': gvcfs_idx, 'doneFiles': doneFiles}

def make_intervals(outputDir, intDir, wildcards, dict_file):
    """Creates interval list files for parallelizing haplotypeCaller and friends. Writes one contig/chromosome per list file."""
    
    with open(dict_file, "r") as f:  # Read dict file to get contig info
        contigs = defaultdict()
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                contigs[chrom] = ln
        
        interval_file = os.path.join(outputDir,wildcards.Organism,wildcards.refGenome,intDir, f"{wildcards.refGenome}_intervals_fb.bed")
        with open(interval_file, "w") as fh:
            for contig, ln in contig.items():
                print(f"{contig}:1-{ln}", file=fh)
        
        for i, (contig, ln) in enumerate(contigs.items()):
            interval_list_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{i}.list")
            with open(interval_list_file, "w") as f:
                print(f"{contig}:1-{ln}", file=f)
            
            