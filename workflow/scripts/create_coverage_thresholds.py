from snakemake.script import snakemake
from snakemake.exceptions import WorkflowError
import math

# read chrom coverage values and compute min/max
cov_thresh = {}
stdv_scale = snakemake.params["cov_threshold_stdev"]
rel_scale = snakemake.params["cov_threshold_rel"]
mean_lower = snakemake.params["cov_threshold_lower"]
mean_upper = snakemake.params["cov_threshold_upper"]

# check that correct settings are set

if stdv_scale:
    if rel_scale:
        raise WorkflowError(
            "Both cov_threshold_stdev and cov_threshold_rel are set, please choose one and make sure the other variable is empty in the config file."
        )
    elif mean_lower or mean_upper:
        raise WorkflowError(
            "Both cov_threshold_stdev and cov_threshold_lower/cov_threshold_upper are set, please choose one and make sure the other variable is empty in the config file."
        )
elif rel_scale:
    if mean_lower or mean_upper:
        raise WorkflowError(
            "Both cov_threshold_rel and cov_threshold_lower/cov_threshold_upper are set, please choose one and make sure the other variable is empty in the config file."
        )
elif mean_lower:
    if not mean_upper:
        mean_upper = 50000
elif mean_upper:
    if not mean_lower:
        mean_lower = 1
else:
    raise WorkflowError(
        "Use coverage filter is True, but you did not specify coverage filtering options in the config. Please check."
    )

with open(snakemake.input["stats"]) as stats:
    for line in stats:
        if "mean" in line:
            continue

        fields = line.split()
        mean = float(fields[1])
        stdev = math.sqrt(mean)
        # 0 is chr, 1 is mean
        if stdv_scale:
            cov_thresh[fields[0]] = {
                "low": mean - (stdev * float(stdv_scale)),
                "high": mean + (stdev * float(stdv_scale)),
            }
        elif rel_scale:
            cov_thresh[fields[0]] = {
                "low": mean / float(rel_scale),
                "high": mean * float(rel_scale),
            }
        else:
            cov_thresh[fields[0]] = {
                "low": float(mean_lower),
                "high": float(mean_upper),
            }

# Write the thresholds to a TSV file
with open(snakemake.output[0], "w") as output_file:
    output_file.write("chrom\tmin\tmax\n")  # Header line, if needed
    for chrom, thresholds in cov_thresh.items():
        if chrom == "total":
            continue
        output_file.write(f"{chrom}\t{thresholds['low']}\t{thresholds['high']}\n")
