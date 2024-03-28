import glob
import re
import sys
import os
from pathlib import Path

# Get utils. This is not great, but we can move to setup.py and install via pip later if want
utils_path = (Path(workflow.main_snakefile).parent.parent.parent).resolve()
if str(utils_path) not in sys.path:
    sys.path.append(str(utils_path))

import pandas as pd
import snparcher_utils

samples = snparcher_utils.parse_sample_sheet(config)

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
                raise WorkflowError(f"Reference gff {ref} does not exist")
            elif ref.rsplit(".", 1)[1] == '.gz':
                raise WorkflowError(f"Reference gff {ref} must be unzipped first.")
        return _refs
    else:
        return []
