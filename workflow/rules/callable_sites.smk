
## RULES ##

## plan is to make a callable sites that filters on both mappability and coverage


rule compute_covstats:
    input:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".covstats.bg.gz"
    run:
        import gzip
        with gzip.open(input.bgz) as f:
            with gzip.open(output.cov, 'wb') as covbed:
                for line in f:
                    fields = line.split()
                    if (len(fields) < 4):
                        continue
                    cov_fields = map(int, fields[3:])
                    covsum = sum(cov_fields)
                    mean = covsum / len(fields[3:])
                    print(fields[0], fields[1], fields[2], covsum, mean, file=covbed)
