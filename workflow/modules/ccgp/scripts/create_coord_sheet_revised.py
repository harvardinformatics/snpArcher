import sys
# sys.path.append(
#     ".."
# )
import pandas as pd
import os
from db import get_mongo_client
import math
#these should be native
#import argparse
#mport re

client = get_mongo_client()
db = client["ccgp_dev"]
collection = db["sample_metadata"]
wd_scripts = os.getcwd()
#print(wd_scripts)
#os.chdir("..")
wd_parent = os.getcwd() #this should be ccgp-reruns/

def mongo_get_coords(project, ref_genome, sample_type, vcf_samples):

    data_list = []
    no_coords = []
    sample_bucket = []
    query_criteria = {"ccgp-project-id": project}
    samps_no_coords = vcf_samples.copy()

    sample_entries_for_project = collection.find(query_criteria)

    for sample in sample_entries_for_project:
        sample_name_exist = False
        minicore_seq_exist = False
        # if sample.get("*sample_name", "Womp_Womp") in vcf_samples:
        #     sample_name_exist = True
        if str(sample["*sample_name"]) in vcf_samples:
            sample_name_exist = True

        if sample.get("minicore_seq_id", "Womp") in vcf_samples:
            minicore_seq_exist = True
        if sample_name_exist and minicore_seq_exist:
            print('IDENTICAL SAMPLE_NAME & MINICORE') 
            print('')
            sample_name = str(sample.get("*sample_name", ""))
            print(sample_name)
            lat = sample.get("lat", "")
            long = sample.get("long", "")

            if lat == '' or lat == 0 or str(lat).lower() == 'nan':
                no_coords.append(sample_name)
            else:
                data_list.append({"Sample Name": sample_name, "Long": float(long), "Lat": float(lat)})
            samps_no_coords.remove(sample_name)
            continue

        if sample_name_exist: #IF THE ENTRY USES SAMPLE_NAMES  
            print('SAMPLE_NAME') 
            print('')
            sample_name = str(sample.get("*sample_name", ""))
            lat = sample.get("lat", "")
            long = sample.get("long", "")
            print(str(lat).lower()=='nan')
            print(str(lat).lower())
            if lat == '' or lat == 0 or str(lat).lower() == 'nan':
                no_coords.append(sample_name)
            else:
                data_list.append({"Sample Name": sample_name, "Long": float(long), "Lat": float(lat)})
            print(sample_name)
            samps_no_coords.remove(sample_name)
        #print(sample_name)
        if minicore_seq_exist:
            print("MINICORE")
            print('')
            sample_name = sample.get("minicore_seq_id", "")
            lat = sample.get("lat", "")
            long = sample.get("long", "")

            if lat == '' or lat == 0 or str(lat).lower() == 'nan':
                no_coords.append(sample_name)
            else:
                data_list.append({"Sample Name": sample_name, "Long": float(long), "Lat": float(lat)})
            print(sample_name)
            samps_no_coords.remove(sample_name)
        # else:
        #     no_coords.append(sample)

    df = pd.DataFrame(data_list)
    df_nocoords = pd.DataFrame(no_coords)
    num_rows = df.shape[0]

    if num_rows == 0:
        print(f"No results found for project {project}. Skipping...")

    print(f"Number of samples considered = {len(vcf_samples)}")
    print(f"Number of samples WITH coords for {project}: {num_rows}")
    print(f"Number of samples WITHOUT coords for {project}: {len(no_coords)}")
    output_file_path = os.path.join(wd_scripts, "results", ref_genome, "CCGP", f"{project}.coords.txt")
    output_file_path_nocoords = os.path.join(wd_scripts, "results", ref_genome, "CCGP", f"{project}.no_coords.txt")
    #print(output_file_path)

    df.to_csv(output_file_path, sep="\t", index=False, header=False)
    # for samp in samps_no_coords:
    #     print(samp)
    

    with open(output_file_path_nocoords, 'w') as file:
        for samp in no_coords:
            file.write(samp + '\n')


    print(f'Got {project} coords.')
    print('')

if __name__ == "__main__":
    ccgp_project_id = snakemake.params["project_id"]
    ref_genome = snakemake.params["ref_genome"]
    sample_type = snakemake.params["sample_id"]
    samps = snakemake.input["samps"]

    vcf_samples = []
    coord_samples_file = samps
    with open(coord_samples_file) as file:
        for sample in file:
            vcf_samples.append(sample.strip())

    mongo_get_coords(ccgp_project_id, ref_genome, sample_type, vcf_samples)

    

