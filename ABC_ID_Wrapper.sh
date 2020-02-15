#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/ABC_id'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SLC_FAM=$H/ABC_REF/Input_files/ABC_Families.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
THREADS=36

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





### phylogeny
#mkdir ABC_phylo
#source ./ABC_ID_SCRIPTS/ABC_phylo.sh