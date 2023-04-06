import pymongo
import pandas as pd
import sys
sys.path.append(
    ".."
)
from utils.db import get_mongo_client

# Finding orphans in the database
# Compare the files in the database to the files in the s3 bucket
# If the file is not in mongo then print it out


client = get_mongo_client()
db = client['ccgp_dev']
collection = db['sample_metadata']


#if you are looking for a specific project id use this query (uncomment the line below) and -p "project_id" in the command line
#docs = collection.find({"ccgp-project-id": config['pid']})

filenames = []
docs = collection.find()
for doc in docs:

    if doc.get(f"read1") is not None:
        if isinstance(doc.get(f"read1"), str):
            filenames.append(doc[f"read1"]).strip()
            #filenames.strip() == filenames
    if doc.get(f"read2") is not None:
        if isinstance(doc.get(f"read2"), str):
            filenames.append(doc[f"read2"]).strip()
            #filenames.strip() == filenames

    if doc.get(f"run1_read1") is not None:
        for i in range(1,10):
            if doc.get(f"run{i}_read1") is not None:
                if isinstance(doc.get(f"run{i}_read1"), str):
                    filenames.append(doc[f"run{i}_read1"].strip())
            if doc.get(f"run{i}_read2") is not None:
                if isinstance(doc.get(f"run{i}_read2"), str):
                    filenames.append(doc[f"run{i}_read2"].strip())
    
    if doc.get(f"files") is not None:
        if isinstance(doc.get(f"files"), list):
            for i in doc[f"files"]:
                filenames.append(i.strip())
                #print(i)
    
with open(r'orphan_reads_result/mongo-filenames-file.txt', 'w') as fp:
    fp.write('\n'.join(filenames))
#print(filenames)
print(len(filenames))


rule all:
    input:
        "orphan_reads_result/orphaned_reads.txt"
        
#Find what is missing in mongo-filenames-file.txt from aws-complete-file.txt
rule compare_files:
    input:
        "orphan_reads_result/aws-complete-file.txt",
        "orphan_reads_result/mongo-filenames-file.txt"
    output:
        "orphan_reads_result/orphaned_reads.txt"
    shell:
        "grep -F -x -v -f orphan_reads_result/mongo-filenames-file.txt orphan_reads_result/aws-complete-file.txt > orphan_reads_result/orphaned_reads.txt"

rule mongo_filenames:
    output:
        "orphan_reads_result/mongo-filenames-file.txt"
    shell:
        "orphan_reads_result/mongo-filenames-file.txt"

rule clean_s3_bucket_query:
    input:
        "orphan_reads_result/aws-complete-file-1.txt"
    output:
        "orphan_reads_result/aws-complete-file.txt"
    shell:
        "awk '{{print $4}}' orphan_reads_result/aws-complete-file-1.txt > orphan_reads_result/aws-complete-file.txt"

rule s3_bucket_query:
    output:
        "orphan_reads_result/aws-complete-file-1.txt"
    shell:
        "aws s3 ls s3://ccgp --endpoint=$endpoint | grep .gz > orphan_reads_result/aws-complete-file-1.txt"
        
   

#create a file that contains all the files in the s3 bucket without the first 2 collums (date and time)
#aws s3 ls s3://ccgp --endpoint=$endpoint | grep .fastq.gz | cut -d " " -f4

#Find what is missing in mongo-filenames-file.txt from aws-complete-file.txt
#sort mongo-filenames-file.txt mongo-filenames-file.txt aws-complete-file.txt | uniq -u > orphans_on.txt

