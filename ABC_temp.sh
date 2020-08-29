#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=80
THREADS=4
cd $H 


###### 4) Search proteomes
for i in ./proteomes/*.faa; do
 ./ABC_ID_SCRIPTS/ABC_search/ABC_search.sh -proteome $i -outdir ./ABC_search -threads $THREADS -minlen 250
done
