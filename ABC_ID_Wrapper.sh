#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=75
THREADS=8

### 2) Set Home working Directory
cd $H 

### 3) Prepare proteomes
mkdir proteomes
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh
  
###### 2) Build database materials
source ./ABC_ID_SCRIPTS/ABC_Model_database_build.sh


##### 3) Run BUSCO for quality assessment
source ./ABC_ID_SCRIPTS/ABC_BUSCO.sh

###### 4) Search proteomes
mkdir -p ABC_search
mkdir -p preliminary_ABC preliminary_ABC/proteomes preliminary_ABC/dicts
for i in ./proteomes/*.faa; do
 ./ABC_ID_SCRIPTS/ABC_search.sh -target $i -hmm_profile ./model_database/HMM_databases/Only_ABCs.hmm -outdir ABC_search -threads $THREADS
done

###### 5)Filter based on number of NBDs
Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R

###### 5) CAFE
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees

#source ./ABC_ID_SCRIPTS/ABC_species_phylogeny.sh
cp ./GENERAL_REFERENCE/CAFE/ultrametric_tree_backup/*.nwk ./CAFE/clean_raxml_trees/

#Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_prep.R  
#source ./ABC_ID_SCRIPTS/ABC_CAFE_run_full.sh
#Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_figures.R
  
###### 6) Produce Figures and Tables
#Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R
  

### 7) Make phylogeny for relevant species 

#source ./ABC_ID_SCRIPTS/ABC_phylo.sh
#Rscript ~/Applications/Custom_Applications/ggtree_clean_phylogeny.R #### NEED TO GET ALL TREES IN ONE FOLDER
#cp ./GENERAL_REFERENCE/phylo_premade/* ./phylo/clean_trees/


