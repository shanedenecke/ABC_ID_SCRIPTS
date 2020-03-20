#!/usr/bin/env bash
H=~/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=0
THREADS=14

cd $H

###### 1) proteome prepare
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh

###### 2) Build database materials
source ./ABC_ID_SCRIPTS/ABC_Model_database_build_join.sh

###### 3) Search proteomes
mkdir ABC_search
mkdir preliminary_ABC
mkdir preliminary_ABC/proteomes
mkdir preliminary_ABC/dicts
for i in ./proteomes/*; do
  source ./ABC_ID_SCRIPTS/ABC_search_join.sh $i
done
  
  ###### 4)Filter based on number of NBDs
  Rscript ./ABC_ID_SCRIPTS/ABC_domain_filter_new.R $QUAL_THRESH
  
  ###### 5) CAFE
  mkdir CAFE
  mkdir ./CAFE/clean_raxml_trees
  #source ./ABC_ID_SCRIPTS/ABC_Ultrametric_tree_generate.sh
  cp ./ABC_REF/ultrametric_tree_backup/*.tre ./CAFE/clean_raxml_trees/
  Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_prep.R
  source ./ABC_ID_SCRIPTS/ABC_CAFE_run_full.sh
  Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_figures.R
  
  ###### 6) Produce Figures and Tables
  Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R
  

### 7) Make phylogeny for relevant species 
mkdir phylo
mkdir ./phylo/clean_trees
source ./ABC_ID_SCRIPTS/ABC_phylo.sh
Rscript ./ABC_ID_SCRIPTS/ABC_phylo_ggtree_clean.R $H/phylo/clean_trees/
#cp ./ABC_REF/phylo_premade/* ./phylo/clean_trees/


