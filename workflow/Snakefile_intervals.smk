import os
from collections import defaultdict
import sys
import yaml
import helperFun
import pandas as pd

configfile: "config.yaml"
res_config = yaml.load(open("resources.yaml"))

# create directory to store GATK lists
"""if not os.path.isdir(intDir + "gatkLists"):
    os.system("mkdir -p " + intDir + "gatkLists")"""

samples = pd.read_table(config["samples"], sep=",").replace(' ', '_', regex=True)
org_ref = set(zip(samples.Organism.tolist(), samples.refGenome.tolist()))  # To get only unique combinations of organism and ref accession.
ORGANISM, REFGENOME = map(list, zip(*org_ref))  # Split above back to indiviudal lists for expand. There probably is a better way?

### workflow ###

rule all:
    input:
        expand(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed", zip, Organism=ORGANISM, refGenome=REFGENOME)

include: "rules/intervals.smk"




