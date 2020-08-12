#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=75
THREADS=8

### 2) Set Home working Directory
cd $H 



###### 4) Search proteomes
mkdir -p ABC_search
for i in ./proteomes/*; do
 ./ABC_ID_SCRIPTS/ABC_search.sh -target $i -hmm_profile ./model_database/HMM_databases/Only_ABCs.hmm -outdir ABC_search -threads $THREADS
done

