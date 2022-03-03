
## RULES ##

## plan is to make a callable sites that filters on both mappability and coverage


rule compute_covstats:
    input:
        bgz = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + ".bg.gz"
    output:
        cov = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}" + "covstats.bg.gz"
    run:
        import gzip
        with gzip.open(input.bgz) as f:
            with gzip.open(output.cov, 'wb') as covbed:
                for line in f:
                    fields = line.split()
                    sum = sum(int(fields[3,]))
                    mean = sum / length(fields[3,])
                    print(fields[0], fields[1], fields[2], sum, mean, file=covbed)
