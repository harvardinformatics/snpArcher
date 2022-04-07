from pyd4 import D4File,D4Builder
import math

#read chrom coverage values and compute min/max
cov_thresh = {}
stdv_scale = float(snakemake.params["cov_threshold"])
with open(snakemake.input["stats"]) as stats:
    for line in stats:
        fields=line.split()
        stdev = math.sqrt(float(fields[1]))
        cov_thresh[fields[0]] = {'low' : float(fields[2]) - (stdev * stdv_scale),
            'high' : float(fields[2]) + (stdev * stdv_scale)}

#read d4 file into python, convert to
covfile = D4File(snakemake.input["d4"])
covmat = covfile.open_all_tracks()

with open(snakemake.output["covbed"], mode='w') as covbed:
    for chrom in covfile.chroms():
        for values in covmat.enumerate_values(chrom[0],0,chrom[1]):
            covs=values[2]
            res1=math.fsum(covs)
            if res1 <= cov_thresh[chrom[0]]['high'] and res1 >= cov_thresh[chrom[0]]['low']:
                print(chrom[0], values[1])
