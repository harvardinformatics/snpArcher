#!/bin/bash
export CLUSTER_NAME="snakemake"
export ZONE="us-central1"
export PROJECT_ID="ccgp-ucsc"

gcloud container clusters create $CLUSTER_NAME \
    --project=${PROJECT_ID} \
    --zone=${ZONE} \
    --machine-type="n1-standard-1" \
    --num-nodes=1 \

gcloud container node-pools create "homestar-runner3" \
    --cluster=${CLUSTER_NAME} \
    --machine-type="n2d-standard-32" \
    --num-nodes=0 \
    --zone=${ZONE} \
    --preemptible \
    --project=${PROJECT_ID} \
    --scopes storage-rw \
    --image-type=UBUNTU \
    --disk-size=2TB \
    --disk-type=pd-balanced \
    --enable-autoscaling \
    --min-nodes=0 \
    --max-nodes=100 \