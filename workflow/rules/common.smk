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

def make_intervals(outputDir, intDir, wildcards, dict_file, num_intervals, min_int_size):
    """Creates intervals by splitting all contigs by N number of splits. num_intervals is the number of interval FILES. min_intsize defines
    the bp size of each interval that a contig can be split into."""
    
    with open(dict_file, "r") as f:  # Read dict file to get contig info
        contigs = defaultdict()
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                contigs[chrom] = ln
        
        ln_sum = sum(contigs.values())
        bp_per_interval = ln_sum // int(num_intervals)
        
        int_file = 0
        running_bp_total = 0
        out = deque()
        
        for chrom, ln in contigs.items():
            if ln <= int(min_int_size):
                out.append(f"{chrom}:1-{ln}")
                running_bp_total += ln
                if running_bp_total >= bp_per_interval:
                    interval_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{int_file}.list")
                    with open(interval_file, "a+") as f:
                        for _ in range(len(out)):
                            line = out.popleft()
                            print(line, file=f)
                    int_file += 1
                    running_bp_total = 0
                    
            else:
                interval_ends = list(range(1, ln, int(min_int_size)))
                interval_lists = [[interval_ends[i-1], interval_ends[i]] for i in range(1, len(interval_ends))]  # Make the split ranges eg 1-100, 100-200 ...
                interval_lists[-1][1] = ln  # Sets last value of last interval list to the size of the contig
                for l in interval_lists:
                    if l[0] != 1:  # Increment first num of each range by 1 if its not 1, so no overlapping
                        l[0] += 1
                    out.append(f"{chrom}:{l[0]}-{l[1]}")
                    running_bp_total += (l[1] - l[0])
                    if running_bp_total >= bp_per_interval:
                        interval_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{int_file}.list")
                        with open(interval_file, "a+") as f:
                            for _ in range(len(out)):
                                line = out.popleft()
                                print(line, file=f)
                        int_file += 1
                        running_bp_total = 0
        if out:
            interval_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{int_file}.list")
            with open(interval_file, "a+") as f:
                for _ in range(len(out)):
                    line = out.popleft()
                    print(line, file=f)

                        
