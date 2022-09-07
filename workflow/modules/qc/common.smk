def get_coords_if_available(wildcards):
    if 'lat' in samples.columns and 'long' in samples.columns:
        return "results/{refGenome}/QC/{prefix}.coords.txt"
    return []

def check_contig_names(fai, touch_file):
    dffai = pd.read_table(fai, sep='\t', header = None)
    fai_result=pd.to_numeric(dffai[0], errors='coerce').notnull().all()
    if fai_result==True:
        print("QC plots not generated because contig names are numeric and plink does not accept numeric contig names")
    elif fai_result==False:
        with open(touch_file, "w") as writer:
            writer.write("contigs are strings")
