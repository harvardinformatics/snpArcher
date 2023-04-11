from distutils.command.upload import upload
from gdrive import CCGPDrive
from pathlib import Path
import pdfkit


def create_readme(project_folder, prefix, project_name, outfile):
    cmd = f"gdrive download {project_folder}"
    link = f"https://drive.google.com/drive/u/0/folders/{project_folder}"
    erik = f"https://erikenbody.github.io/ccgp/{prefix + '_qc.html'}"

    html = f"""
    <img src="https://raw.githubusercontent.com/cademirch/ccgp_workflow/main/docs/CCGPhorizontalblue.jpeg" width=auto height=100>
    <p>As part of the planned CCGP workflow, we've finished generating genotype data from your whole-genome resequencing data on your project: {project_name}  Below we include the information needed to access the data. As a first step, we suggest you view the quality control dashboard for your variant data:
    {erik}:</p>
    <p>The idea is to use these outputs in the QC dashboard to identify if there are any clear outliers in your data and to help direct your filtering and downstream analysis plans. Keep in mind these plots are not final analyses, they are just there to give you a sense of the data. In the accompanying email with the data delivery, our bioinformatics team describes anything they noticed in particular with your dataset. Please see the figure descriptions for a general idea of how to use them to identify problematic samples.</p>
    <p>We have hosted your data on Google Drive in:
    {link}</p>
    <p>This folder structure looks like this:</p>
    <pre><code>Project-folder
        ├── QC-Outputs
        ├── Project.filtered.vcf.gz
        └── Project.final.vcf.gz
    </code></pre>
    <p>You can find two VCFs in the top level folder, as well as all the data that was used to generate the QC plots. The small VCF file in the QC folder was used to generate the QC plots, and this dataset has basic quality filters applied and has been coarsely pruned to eliminate linked variants. This file may be a good starting point for exploring your data.</p>
    <p>The &quot;final&quot; VCF corresponds to the raw output of our variant calling pipeline with some basic GATK filters annotated. The &quot;filtered&quot; VCF has removed variants not passing that filter and also samples less than 2x sequencing depth.</p>
    <p>To download the outputs directly to your server we recommend using <a href="https://github.com/prasmussen/gdrive">gDrive</a>. Below you will find the command to download your project results folder. Alternatively, you can download the files to your local machine and then move them to your server using a tool such as rsync or scp.</p>
    <p><code>{cmd}</code></p>
    <p>For a detailed description of the pipeline, please refer to our GitHub README and associated code:
    <a href="https://github.com/cademirch/ccgp_workflow">https://github.com/cademirch/ccgp_workflow</a></p>
    <p>Please do not hesitate to contact members of the bioinformatics team if you have any questions:
    <br>Erik Enbody (<a href="mailto:erik.enbody@gmail.com">erik.enbody@gmail.com</a>)
    <br>Mara Baylis (<a href="mailto:mbaylis@ucsc.edu">mbaylis@ucsc.edu</a>)
    <br>Cade Mirchandani (<a href="mailto:cmirchan@ucsc.edu">cmirchan@ucsc.edu</a>)</p>
    """

    html_header = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
    </head>
    """

    html = html_header + html + "</html>"

    pdfkit.from_string(html, outfile)
    return Path(outfile)


def setup_dirs(prefix: str, drive_client: CCGPDrive) -> dict:
    """Sets up directories in google drive for project's results.

    Args:
        prefix (str): The prefix for the project.
        drive_client (CCGPDrive): CCGPDrive class to interact with Drive.

    Returns:
        dict:
            key: folder
            value: drive id
    """
    # check if folders exist
    current_proj_folders = drive_client.list_files_from_folder("Project Results")
    if prefix not in [f["name"] for f in current_proj_folders]:
        print(
            "Did not find exisitng results folder for this project in Drive. Creating one..."
        )
        project_folder = drive_client.create_folder(
            prefix, drive_client.get_folder_id("Project Results")
        )
        qc_folder = drive_client.create_folder("QC Outputs", project_folder)
        gvcf_folder = drive_client.create_folder("GVCFs", project_folder)
        #bam_folder = drive_client.create_folder("BAMS", project_folder)
        return {
            "project_folder": project_folder,
            "qc_folder": qc_folder,
            "gvcf_folder": gvcf_folder,
        }
        #"bam_folder": bam_folder, (not using bam files)
    else:  # a project results folder for this project DOES exist
        for f in current_proj_folders:
            if f["name"] == prefix:
                project_folder = f["id"]
                print(f"Found folder: {project_folder} for project: {prefix}")
        current_files = drive_client.list_files_from_folder(project_folder, _id=True)

        for f in current_files:
            if f["name"] == "QC Outputs":
                qc_folder = f["id"]
            if f["name"] == "GVCFs":
                gvcf_folder = f["id"]
            # if f["name"] == "BAMS":
                # bam_folder = f["id"]

        if "QC Outputs" not in [f["name"] for f in current_files]:
            print("Didn't find QC folder in exisiting results folder, creating one...")
            qc_folder = drive_client.create_folder("QC Outputs", project_folder)
        if "GVCFs" not in [f["name"] for f in current_files]:
            print(
                "Didn't find GVCF folder in exisiting results folder, creating one..."
            )
            gvcf_folder = drive_client.create_folder("GVCFs", project_folder)
        if "BAMS" not in [f["name"] for f in current_files]:
            print(
                "Didn't find BAMS folder in exisiting results folder, creating one..."
            )
            # bam_folder = drive_client.create_folder("BAMS", project_folder)
        return {
            "project_folder": project_folder,
            "qc_folder": qc_folder,
            "gvcf_folder": gvcf_folder,
            # "bam_folder": bam_folder,
        }


def main():
    drive = CCGPDrive()
    prefix = snakemake.config["final_prefix"]
    dirs = setup_dirs(prefix=prefix, drive_client=drive)
    if snakemake.rule == "readme":
        readme = create_readme(
            project_folder=dirs["project_folder"],
            prefix=prefix,
            project_name=prefix,
            outfile=snakemake.output[0],
        )
        drive.upload_file(readme, folder_id=dirs["project_folder"])
    if snakemake.rule == "upload_qc":
        files_to_upload = [
            "eigenvec",
            "eigenval",
            "depth",
            "dist",
            "distid",
            "king",
            "miss",
            "admix3",
            "admix2",
            "snpqc",
            "faiResult",
            "bed",
            "bim",
            "fam",
            "sumstats",
            "summ",
            "het",
            "fai",
            "qcpdf",
        ]
        for k in files_to_upload:
            drive.upload_file(Path(snakemake.input[k]), folder_id=dirs["qc_folder"])

    if snakemake.rule == "upload_gvcfs":
        drive.upload_file(Path(snakemake.input["gvcf"]), folder_id=dirs["gvcf_folder"])
        drive.upload_file(
            Path(snakemake.input["gvcf_idx"]), folder_id=dirs["gvcf_folder"]
        )
    #if snakemake.rule == "upload_bams":
        #drive.upload_file(Path(snakemake.input["bam"]), folder_id=dirs["bam_folder"])
        #drive.upload_file(Path(snakemake.input["bai"]), folder_id=dirs["bam_folder"])

    if snakemake.rule == "upload_vcfs":
        drive.upload_file(
            Path(snakemake.input["vcf"]), folder_id=dirs["project_folder"]
        )
        drive.upload_file(
            Path(snakemake.input["tbi"]), folder_id=dirs["project_folder"]
        )
        drive.upload_file(
            Path(snakemake.input["filtered"]), folder_id=dirs["project_folder"]
        )
        drive.upload_file(
            Path(snakemake.input["filtered_idx"]), folder_id=dirs["project_folder"]
        )


if __name__ == "__main__":
    main()
