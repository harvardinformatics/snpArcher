import glob
import re
import os
import sys
import pandas as pd
from snakemake.exceptions import WorkflowError

samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)

def get_ref(wildcards):
    if 'refPath' in samples.columns:
        _refs = samples.loc[(samples['refGenome'] == wildcards.refGenome)]['refPath'].dropna().unique().tolist()
        for ref in _refs:
            if not os.path.exists(ref):
                raise WorkflowError(f"Reference genome {ref} does not exist")
            elif ref.rsplit(".", 1)[1] == '.gz':
                raise WorkflowError(f"Reference genome {ref} must be unzipped first.")
        return _refs
    else:
        return []

def get_gff(wildcards):
    if 'refGFF' in samples.columns:
        _refs = samples.loc[(samples['refGenome'] == wildcards.refGenome)]['refGFF'].dropna().unique().tolist()
        for ref in _refs:
            if not os.path.exists(ref):
                raise WorkflowError(f"Reference genome {ref} does not exist")
            elif ref.rsplit(".", 1)[1] == '.gz':
                raise WorkflowError(f"Reference genome {ref} must be unzipped first.")
        return _refs
    else:
        return []
