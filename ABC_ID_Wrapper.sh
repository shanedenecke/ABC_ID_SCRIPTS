#!/usr/bin/env bash
H=/mnt/disk/shane/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.3
THREADS=10

cd $H 

###### 1) proteome prepare
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh
#source ./ABC_ID_SCRIPTS/ABC_proteome_add.sh $H/ABC_REF/non_model_proteomes/append/*

###### 2) Build database materials
source ./ABC_ID_SCRIPTS/ABC_Model_database_build_join.sh


##### 2.5) Run BUSCO for quality assessment
source ./ABC_ID_SCRIPTS/ABC_BUSCO.sh

###### 3) Search proteomes
mkdir ABC_search
mkdir preliminary_ABC
mkdir preliminary_ABC/proteomes
mkdir preliminary_ABC/dicts
for i in ./proteomes/*; do
  source ./ABC_ID_SCRIPTS/ABC_search_join.sh $i
done

###### 4)Filter based on number of NBDs
Rscript ./ABC_ID_SCRIPTS/ABC_domain_filter_new.R

###### 5) CAFE
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees

#source ./ABC_ID_SCRIPTS/ABC_species_phylogeny.sh
cp ./ABC_REF/ultrametric_tree_backup/*.nwk ./CAFE/clean_raxml_trees/

Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_prep.R  
source ./ABC_ID_SCRIPTS/ABC_CAFE_run_full.sh
Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_figures.R
  
###### 6) Produce Figures and Tables
Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R
  

### 7) Make phylogeny for relevant species 

source ./ABC_ID_SCRIPTS/ABC_phylo.sh
#Rscript ~/Applications/Custom_Applications/ggtree_clean_phylogeny.R #### NEED TO GET ALL TREES IN ONE FOLDER
#cp ./ABC_REF/phylo_premade/* ./phylo/clean_trees/


