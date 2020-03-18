#!/usr/bin/env bash
H=~/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.3
THREADS=14

cd $H

###############1) proteome prepare
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh

##############2) Build database materials
source ./ABC_ID_SCRIPTS/ABC_Model_database_build_join.sh

##############3) Search proteomes
mkdir ABC_search
mkdir preliminary_ABC
mkdir preliminary_ABC/proteomes
mkdir preliminary_ABC/dicts
for i in ./proteomes/*; do
  source ./ABC_ID_SCRIPTS/ABC_search_join.sh $i
done



###### Filter
source ./ABC_ID_SCRIPTS/ABC_domain_filter.sh

#### CAFE
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees

#source ./ABC_ID_SCRIPTS/ABC_Ultrametric_tree_generate.sh
cp ./ABC_REF/ultrametric_tree_backup/*.tre ./CAFE/clean_raxml_trees/
Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_prep.R
source ./ABC_ID_SCRIPTS/ABC_CAFE_run_full.sh
Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_figures.R


### phylogeny
mkdir phylo
mkdir ./phylo/clean_trees
source ./ABC_ID_SCRIPTS/ABC_phylo.sh ###################### NEED TO EDIT
#cp ./ABC_REF/phylo_premade/* ./phylo/clean_trees/




