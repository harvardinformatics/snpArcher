import os
from utils.db import get_mongo_client
import pymongo
import datetime
import math
from pprint import pprint
import pandas as pd
import re
import copy
from thefuzz import process
from pathlib import Path
from datetime import datetime
import argparse
from itertools import chain
from utils.gsheets import WGSTracking
from utils.gdrive import CCGPDrive

pd.options.mode.chained_assignment = None


def preprocess_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Processes dataframe from MongoDB to be turned into workflow sheet or SRA sheet."""
    list_of_dicts = []
    for _, row in df.iterrows():
        record = row.to_dict()
        # print(record.get("*sample_name"))
        # print(record.get("files"))
        if record.get(f"run1_read1") is not None:
            pairs = {f"run{i}": None for i in range(1,107)}
            for run in pairs.keys():
                f1, f2 = record.get(f"{run}_read1"), record.get(f"{run}_read2")
                if f1 is None and f2 is None:
                    break
                if not isinstance(f1, str) and not isinstance(f2, str):
                    if math.isnan(f1) and math.isnan(f2):
                        print(f"Found nan file for sample {record.get('*sample_name')}")
                        break
                
                if (f1 is None and f2 is not None) or (f1 is not None and f2 is None):
                    raise ValueError(f"Sample {record['*sample_name']} has missing read pair ({run}): fq1: {f1}, fq2: {f2}")
                if not isinstance(f1, str) or not isinstance(f2, str):
                    raise ValueError(f"Sample {record['*sample_name']} has non string read names ({run}): fq1: {f1}, fq2: {f2}")
                if not f1 or not f2:
                    raise ValueError(f"Sample {record['*sample_name']} has missing read pair ({run}): fq1: {f1}, fq2: {f2}")
                f1, f2 = f1.strip(), f2.strip()
                record_copy = copy.deepcopy(record)
                record_copy["fq1"] = f1
                record_copy["fq2"] = f2
                record_copy["library_ID"] = str(record_copy["*sample_name"]) + f"_{run}"
                record_copy[
                    "title"
                ] = f"Whole genome sequencing of {record_copy['*organism']}"
                list_of_dicts.append(record_copy)
            # continue
                
        if record.get(f"filename1") is not None:
            pairs = {}
            for i in range(1,107,2):
                pairs[i] = [record.get(f"filename{i}"), record.get(f"filename{i+1}")]
            for key, pair in pairs.items():
                f1, f2 = pair
                if f1 is None or f2 is None:
                    continue
                elif not isinstance(f1, str) or not isinstance(f2, str):
                    continue
                elif not f1 or not f2:
                    continue
                else:
                    record_copy = copy.deepcopy(record)
                    record_copy["fq1"] = f1
                    record_copy["fq2"] = f2
                    record_copy["library_ID"] = str(record_copy["*sample_name"]) + f"_{key}"
                    record_copy[
                        "title"
                    ] = f"Whole genome sequencing of {record_copy['*organism']}"
                    list_of_dicts.append(record_copy)
            continue
        '''
        if record.get("files") is not None:
            files = record["files"]
            if isinstance(record.get(f"files"), list):
                for i in record[f"files"]:
                    record.append(i.strip())
                    print(i)
        '''
                
        if record.get("files") is not None:
            
            files = record["files"]
            if isinstance(record.get(f"files"), list):
                files = record[f"files"]
            
                if len(files) == 2:

                    record["fq1"] = files[0]
                    record["fq2"] = files[1]
                    record["library_ID"] = record["fq1"].split("_R1")[0]
                    record["title"] = f"Whole genome sequencing of {record['*organism']}"
                    list_of_dicts.append(record)
            
                elif len(files) >= 4:

                    prefixes = set()
                    for file in files:
                        prefixes.add(
                            re.split("_R\d[._]", file)[0]
                        )  # this splits on the _R1/R2 of the filename

                # Using fuzzy matching here for the following example case:
                # sample has two sets of reads: samp_a1_R1.fq.gz and samp_a1_L001_R1.fq.gz
                # Splitting on R1 gives two prefixes, samp_a1_R1 and samp_a1_L001, the former is a substring of the latter,
                # So something like [f for f in files if p in f for p in prefixes] doesn't work to find pairs of reads.
                    pairs = []
                    for p in prefixes:
                        matches = process.extract(p, files, limit=2)
                        pairs.append([matches[0][0], matches[1][0]])

                    for pair in pairs:
                        record_copy = copy.deepcopy(record)
                        record_copy["fq1"] = pair[0]
                        record_copy["fq2"] = pair[1]
                        record_copy["library_ID"] = record_copy["fq1"].split("_R1")[0]
                        record_copy[
                            "title"
                        ] = f"Whole genome sequencing of {record_copy['*organism']}"
                        list_of_dicts.append(record_copy)
            
            
    out_df = pd.DataFrame(list_of_dicts)
    return out_df


def create_workflow_sheet(
    project_id: str, db_client: pymongo.MongoClient, out_path=None
) -> None:
    """Creates workflow csv for given project_id. If project does not have ref genome accession, a placeholder will be printed."""
    db = db_client["ccgp_dev"]
    ccgp_samples = db["sample_metadata"]
    ccgp_workflow_progress = db["workflow_progress"]
    df = preprocess_dataframe(
        pd.DataFrame(list(ccgp_samples.find({"ccgp-project-id": project_id})))
    )
    num_samples = len(list(ccgp_samples.find({"ccgp-project-id": project_id})))
    print(f"Created sheet for {num_samples} from project: {project_id}")
    
    
    try:
        workflow_df = df[
            [
                "*sample_name",
                "library_ID",
                "ref_genome_accession",
                "*organism",
                "fq1",
                "fq2",
                "lat",
                "long",
                "ccgp-project-id",
            ]
        ]
    except KeyError:
        workflow_df = df[
            [
                "*sample_name",
                "library_ID",
                "*organism",
                "fq1",
                "fq2",
                "lat",
                "long",
                "ccgp-project-id",
            ]
        ]
        workflow_df["ref_genome_accession"] = "refGenomePlaceholder"

    workflow_df["ref_genome_accession"] = workflow_df["ref_genome_accession"].replace(
        "NaN", "refGenomePlaceholder"
    )

    workflow_df["lat"] = workflow_df["lat"].astype(str).str.replace('"', "")
    workflow_df["long"] = workflow_df["long"].astype(str).str.replace('"', "")
    workflow_df["lat"] = workflow_df["lat"].astype(str).str.strip()
    workflow_df["long"] = workflow_df["long"].astype(str).str.strip()
    rename_col_dict = {
        "*sample_name": "BioSample",
        "library_ID": "LibraryName",
        "ref_genome_accession": "refGenome",
        "ccgp-project-id": "Organism",  # using project-id here to coerce organism to be the same for all samps in a project.
        "lat": "lat",
        "long": "long",
    }
    workflow_df = workflow_df.rename(columns=rename_col_dict).drop(
        columns=["*organism"]
    )
    workflow_df["Run"] = workflow_df["LibraryName"]
    workflow_df["BioProject"] = workflow_df["Organism"]
    workflow_df = workflow_df.drop_duplicates()
    filename = f"{project_id}_workflow.csv"
    if out_path is not None:
        file_path = Path(out_path, filename)
    else:
        file_path = Path(filename)
    workflow_df.to_csv(file_path, index=False)
    ccgp_workflow_progress.update_one(
        filter={"project_id": project_id},
        update={"$set": {"workflow_sheet_created": datetime.utcnow()}},
        upsert=True,
    )
    return file_path


def create_sra_sheet(project_id: str, db_client: pymongo.MongoClient) -> None:
    db = db_client["ccgp_dev"]
    ccgp_samples = db["sample_metadata"]
    ccgp_workflow_progress = db["workflow_progress"]
    sra_df = preprocess_dataframe(
        pd.DataFrame(list(ccgp_samples.find({"ccgp-project-id": project_id})))
    )
    sra_df["library_strategy"] = "WGS"
    sra_df["library_source"] = "GENOMIC"
    sra_df["library_selection"] = "RANDOM"
    sra_df["library_layout"] = "PAIRED"
    sra_df["platform"] = "Illumina"
    sra_df["instrument_model"] = "Illumina NovaSeq 6000"
    sra_df["filetype"] = "fastq"
    sra_df = sra_df.rename(columns={"library_prep_method": "design_description"})
    sra_df = sra_df.rename(columns={"*sample_name": "sample_name"})
    sra_cols = [
        "sample_name",
        "library_ID",
        "title",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "platform",
        "instrument_model",
        "design_description",
        "filetype",
        "filename",
        "filename2",
    ]
    sra_df = sra_df[sra_cols]
    filename = f"{project_id}_sra.tsv"
    filepath = Path("sra_sheets", filename)
    sra_df.to_csv(filepath, index=False, sep="\t")
    drive = CCGPDrive()
    drive.upload_file(filepath, folder_name="SRA Submission Sheets")
    ccgp_workflow_progress.update_one(
        filter={"project_id": project_id},
        update={"$set": {"sra_sheet_created": datetime.utcnow()}},
        upsert=True,
    )


def create_biosample_sheet(project_id: str, db_client: pymongo.MongoClient) -> None:

    wgs_sheet = WGSTracking()

    taxon_group = wgs_sheet.project_type()[project_id]

    db = db_client["ccgp_dev"]
    ccgp_samples = db["sample_metadata"]
    ccgp_workflow_progress = db["workflow_progress"]
    df = pd.DataFrame(list(ccgp_samples.find({"ccgp-project-id": project_id})))
    df = df.rename(columns={"library_prep_method": "design_description"})

    if "Preferred Sequence ID" in df.columns:
        df = df.rename(columns={"Preferred Sequence ID": "sample_title"})

    df["lat_lon"] = df["lat"].astype(str) + "," + df["long"].astype(str)

    df["isolate"] = (
        df["*organism"].astype(str).replace(" ", "_", regex=True)
        + "_"
        + df["*sample_name"]
    )  # to avoid the issue where ncbi complains about non-unique records.
    df["bioproject_accession"] = ""
    plant_cols = [
        "*sample_name",
        "sample_title",
        "bioproject_accession",
        "*organism",
        "isolate",
        "cultivar",
        "ecotype",
        "age",
        "*gen_loc_name",
        "*tissue",
        "biomaterial_provider",
        "cell_line",
        "cell_type",
        "collected_by",
        "collection_date",
        "culture_collection",
        "disease",
        "disease_stage",
        "genotype",
        "growth_protocol",
        "isolation_source",
        "lat_lon",
        "phenotype",
        "population",
        "sample_type",
        "sex",
        "specimen_voucher",
        "temp",
        "treatment",
        "description",
        "library_prep_method",
    ]

    invert_cols = [
        "*sample_name",
        "sample_title",
        "bioproject_accession",
        "*organism",
        "isolate",
        "breed",
        "host",
        "isolation_source",
        "*collection_date",
        "*geo_loc_name",
        "*tissue",
        "age",
        "altitude",
        "biomaterial_provider",
        "collected_by",
        "depth",
        "dev_stage",
        "env_broad_scale",
        "host_tissue_sampled",
        "identified_by",
        "lat_lon",
        "sex",
        "specimen_voucher",
        "temp",
        "description",
        "library_prep_method",
    ]

    animals_cols = [
        "*sample_name",
        "sample_title",
        "bioproject_accession",
        "*organism",
        "strain",
        "isolate",
        "breed",
        "cultivar",
        "ecotype",
        "age",
        "dev_stage",
        "*sex",
        "*tissue",
        "biomaterial_provider",
        "birth_date",
        "birth_location",
        "breeding_history",
        "breeding_method",
        "cell_line",
        "cell_subtype",
        "cell_type",
        "collected_by",
        "collection_date",
        "culture_collection",
        "death_date",
        "disease",
        "disease_stage",
        "genotype",
        "geo_loc_name",
        "growth_protocol",
        "health_state",
        "isolation_source",
        "lat_lon",
        "phenotype",
        "sample_type",
        "specimen_voucher",
        "store_cond",
        "stud_book_number",
        "treatment",
        "description",
        "library_prep_method",
    ]
    file_path = Path("biosample_sheets", f"{project_id}_biosample.tsv")
    if taxon_group == "Vertebrate":
        df = df[df.columns.intersection(animals_cols)]
        df.to_csv(file_path, index=False, sep="\t")
    elif taxon_group == "Invertebrate":
        df = df[df.columns.intersection(invert_cols)]
        df.to_csv(file_path, index=False, sep="\t")
    elif taxon_group == "Plant":
        df = df[df.columns.intersection(plant_cols)]
        df.to_csv(file_path, index=False, sep="\t")
    else:
        raise ValueError(f"Unexpected taxon group: {taxon_group}")

    drive = CCGPDrive()
    drive.upload_file(file_path, folder_name="BioSample Submission Sheets")


def main():
    db = get_mongo_client()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        dest="project_id",
        required=True,
        help="Project id to generate sheets for.",
    )
    subparser = parser.add_subparsers(
        dest="command", title="subcommands", description="valid subcommands"
    )
    metadata = subparser.add_parser("workflow", description="Create workflow sheet")
    sra = subparser.add_parser("sra", description="Create SRA submission sheet")
    biosample = subparser.add_parser(
        "biosample", description="Create BioSample submission sheet"
    )
    both = subparser.add_parser(
        "all", description="Create all sheets (worfklow, sra, biosample)"
    )
    args = parser.parse_args()

    if args.command == "workflow":
        create_workflow_sheet(args.project_id, db)
    elif args.command == "sra":
        create_sra_sheet(args.project_id, db)
    elif args.command == "biosample":
        create_biosample_sheet(args.project_id, db)
    elif args.command == "all":
        create_workflow_sheet(args.project_id, db)
        create_sra_sheet(args.project_id, db)
        create_biosample_sheet(args.project_id, db)


if __name__ == "__main__":
    main()
