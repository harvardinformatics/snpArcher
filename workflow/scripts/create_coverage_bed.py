from pyd4 import D4File,D4Builder
import math

#read chrom coverage values and compute min/max
cov_thresh = {}
stdv_scale = snakemake.params["cov_threshold_stdev"]
rel_scale = snakemake.params["cov_threshold_rel"]
mean_lower = snakemake.params["cov_threshold_lower"]
mean_upper = snakemake.params["cov_threshold_upper"]

#check that correct settings are set

if stdv_scale:
    if rel_scale:
        raise WorkflowError(f"Both cov_threshold_stdev and cov_threshold_rel are set, please choose one and make sure the other variable is empty in the config file.")
    elif mean_lower or mean_upper:
        raise WorkflowError(f"Both cov_threshold_stdev and cov_threshold_lower/cov_threshold_upper are set, please choose one and make sure the other variable is empty in the config file.")
elif rel_scale:
    if mean_lower or mean_upper:
         raise WorkflowError(f"Both cov_threshold_rel and cov_threshold_lower/cov_threshold_upper are set, please choose one and make sure the other variable is empty in the config file.")
elif mean_lower:
    if not mean_upper:
        mean_upper = 50000
elif mean_upper:
    if not mean_lower:
        mean_lower = 1
else:
    raise WorkflowError(f"Use coverage filter is True, but you did not specify coverage filtering options in the config. Please check.")

with open(snakemake.input["stats"]) as stats:
    for line in stats:
        if "mean" in line:
            continue

        fields=line.split()
        mean = float(fields[1])
        stdev = math.sqrt(mean)
        #0 is chr, 1 is mean
        if stdv_scale:
            cov_thresh[fields[0]] = {
                'low' : mean - (stdev * float(stdv_scale)),
                'high' : mean + (stdev * float(stdv_scale))
                }
        elif rel_scale: 
            cov_thresh[fields[0]] = {
                'low' : mean / float(rel_scale),
                'high' : mean * float(rel_scale)
                }
        else:
            cov_thresh[fields[0]] = {
                'low' : float(mean_lower),
                'high' : float(mean_upper)
                }

#read d4 file into python, convert to
covfile = D4File(snakemake.input["d4"])
covmat = covfile.open_all_tracks()

with open(snakemake.output["covbed"], mode='w') as covbed:
    for chrom in covfile.chroms():
        for values in covmat.enumerate_values(chrom[0],0,chrom[1]):
            covs=values[2]
            #get mean coverage for window
            res1=math.fsum(covs)/len(covs)
            if res1 <= cov_thresh[chrom[0]]['high'] and res1 >= cov_thresh[chrom[0]]['low']:
                print(chrom[0], values[1], values[1]+1, file=covbed, sep="\t")
