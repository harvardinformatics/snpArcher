from os import sched_setscheduler
from utils.create_sheets import create_workflow_sheet
from utils.db import get_mongo_client
from google.cloud import storage
from os.path import join
import pymongo
from git import Repo
import pathlib
from pathlib import Path
import snakemake
import argparse
from snakemake.resources import DefaultResources
import yaml
import os


def check_for_reads(project_id: str, db_client: pymongo.MongoClient) -> None:
    db = db_client["ccgp_dev"]
    ccgp_samples = db["sample_metadata"]
    
    docs = list(
        ccgp_samples.find(
            {"ccgp-project-id": project_id}, {"*sample_name": 1, "files": 1}
        )
    )  # get all docs for given project-id, and only get files and sample field
    

    filenames = []
    for d in docs:
        try:
            filenames.extend(d["files"])
        except KeyError:  # sample doesn't have files
            print(
                f"{d['*sample_name']} does not have files in the database, check that this is correct."
            )
    storage_client = storage.Client()
    bucket = storage_client.bucket("ccgp-raw-reads")
    for f in filenames:
        f = join(project_id, f)
        file = bucket.blob(f)
        if not file.exists():
            print(f"WARNING: File: '{f}' not found in bucket.")


def create_sheet(project_id: str, db_client) -> Path:
    """
    This is where I (mara) am making edits 10/5/22
    Need to prevent overwriting


    Creates workflow sheet for project. Will overwrite existing sheet without warning.

    Args:
        project_id (str): The project id
        db_client (_type_): Pymongo client

    Raises:
        ValueError: If 'refGenomePlaceholder' is in the workflow sheet error is raised. With sed command to fix.

    Returns:
        Path: Path to sample sheet created.
    """
    repo = Repo()  # initialize Repo class

    path = Path(f"sample_sheets/{project_id}_workflow.csv")
    if path.exists():
        print(
            f"Sample sheet: '{path}' already exists! Not making a new one. Delete this one and rerun if you want a new one."
        )
        with open(path, "r") as f:
            for line in f:
                if "refGenomePlaceholder" in line:
                    print(
                        f"'refGenomePlaceholder' found in existing workflow sheet: '{path}'. You can correct this by running the following command:\n sed -i 's|refGenomePlaceholder|<accession>|g' f{path}. Make sure to 'git add -f {path}' after the sed command."
                    )
                    break
        
        print(f"Make sure to run 'git add {str(path)}'")
        return path
    sheet = create_workflow_sheet(project_id, db_client, "sample_sheets")
    
    
    # check if 'refGenomePlaceholder' in sheet
    with open(sheet, "r") as f:
        for line in f:
            if "refGenomePlaceholder" in line:
                print(
                    f"'refGenomePlaceholder' found in workflow sheet: '{sheet}'. You can correct this by running the following command:\n sed -i 's|refGenomePlaceholder|<accession>|g' f{sheet}."
                )
                break
    
    print(f"Make sure to run 'git add {str(sheet)}'")
    
    return sheet


def create_snakemake_config(project_id: str, sheet: Path):
    repo = Repo()
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    d = {
        "samples": str(sheet),
        "resource_config": "config/resources.yaml",
        "final_prefix": project_id,
        "intervals": False,
        "sentieon": True,
        "sentieon_lic": "10.128.0.63:8990",
        "remote_reads": True,
        "remote_reads_prefix": f"ccgp-raw-reads/{project_id}",
        "cov_filter":True,
        "CCGP":False,
        "cov_threshold_stdev": 2,
        #"cov_threshold_lower": ,
        #"cov_threshold_upper": ,
    }
    for k in d.keys():
        config[k] = d[k]
    # print(config)
    path = Path(f"config/{project_id}_config.yaml")
    if not path.exists():
        with open(path, "w") as out:
            out.write(yaml.dump(config, default_flow_style=False))
        # repo.index.add(str(path))
    else:
        print(
            f"Config: '{path}' already exists. Not overwriting. Adding it to git though, just in case."
        )
        # repo.index.add(str(path))
    print(f"Make sure to run 'git add {str(path)}'")
    return d


def create_snakemake_opts(
    project_id,
    config: dict,
    dry_run=False,
    qc_only=False,
    upload_only=False,
) -> dict:

    if qc_only:
        snakefile = "workflow/modules/qc/Snakefile"
        profile = "profiles/gls-sentieon"
    if upload_only:
        snakefile = "workflow/upload.smk"
        profile = "profiles/upload"
    else:
        snakefile = "workflow/Snakefile"
        profile = "profiles/gls-sentieon"

    d = {
        "snakefile": snakefile,
        "config": f"config_file=config/{project_id}_config.yaml",
        "default-remote-prefix": f"ccgp-workflow-results/{project_id}",
        "profile": profile,
    }
    return d


def copy_license(project_id):
    storage_client = storage.Client()
    bucket = storage_client.bucket("ccgp-workflow-results")
    og_license_file = bucket.blob("UCSC_Enbody_eval.lic")
    copy_license_file = bucket.blob(f"{project_id}/UCSC_Enbody_eval.lic")
    if not copy_license_file.exists():
        copy = bucket.copy_blob(
            blob=og_license_file,
            destination_bucket=bucket,
            new_name=copy_license_file.name,
        )


def main():
    db_client = get_mongo_client()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        dest="project_id",
        required=True,
        help="Project id.",
    )
    parser.add_argument(
        "-n", dest="dry_run", action="store_true", help="Do a dry run, duh."
    )
    parser.add_argument(
        "--qc-only",
        dest="qc_only",
        action="store_true",
        help="Only run the qc workflow.",
    )
    parser.add_argument(
        "--upload-only",
        dest="upload_only",
        action="store_true",
        help="Only run the upload workflow.",
    )

    args = parser.parse_args()

    if args.upload_only and args.qc_only:
        raise ValueError(
            "Cannot supply '--qc_only' and '--upload_only'. Please pick one."
        )
    if not args.upload_only or not args.qc_only:
        check_for_reads(args.project_id, db_client)
    sheet = create_sheet(args.project_id, db_client)
    config = create_snakemake_config(args.project_id, sheet)
    opts = create_snakemake_opts(
        args.project_id, config, args.dry_run, args.qc_only, args.upload_only
    )
    #copy_license(args.project_id) #no longer needed

    o = [f"--{k} {v}" for k, v in opts.items()]
    print("snakemake", *o, sep=" ")


if __name__ == "__main__":
    main()
