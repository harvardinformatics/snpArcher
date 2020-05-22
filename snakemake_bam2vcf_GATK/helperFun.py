#!/usr/bin/python -tt

def makeMapFilesForGenomicsDBImport(SAMPLES, LISTS, dbDir, gvcfDir):
	for l in LISTS:
		f=open(dbDir + "DB_mapfile" + l, 'w')
		for s in SAMPLES:
			print(s, gvcfDir + s+"_L"+l+".raw.g.vcf", sep="\t", file=f)
		f.close()

