from snakemake.exceptions import WorkflowError

"""
Reads output file from ScatterIntervalsByNs and puts intervals into (relatively) specified num of equal groups.
Writes intervals groups to individual files for use by HaplotypeCaller
"""


def make_intervals(in_file: str, num_intervals: int, output_files: list[str]) -> None:

    intervals = []

    with open(in_file, "r") as f:
        for line in f:
            if not line.startswith("@"):
                line = line.strip().split()
                chrom, start, end, = line[0], int(line[1]), int(line[2])
                size = end - start
                intervals.append((chrom, start, end, size))

    if num_intervals > len(intervals):
        if in_file[:2] == 'db':
            raise WorkflowError("Number of intervals for db creation is greater than actual number of intervals created by Picard.")
        else:
            raise WorkflowError("Number of intervals for gvcf creation is greater than actual number of intervals created by Picard.")

    groups = [[] for i in range(num_intervals)]
    sums = {i: 0 for i in range(num_intervals)}
    c = 0
    for chrom, start, end, size in sorted(intervals, key=lambda x: x[3]):
        for i in sums:
            if c == sums[i]:
                groups.append((chrom, start, end))
                break
        sums[i] += size
        c = min(sums.values())

    for file in output_files:
        with open(file, "w") as f:
            chrom, start, end = groups.pop()
            print(f"{chrom}:{start}-{end}", file=f)

def main():
    make_intervals(snakemake.input["in_file"], snakemake.params["max_intervals"], snakemake.output["out_files"])

if __name__ == "__main__":
    main()