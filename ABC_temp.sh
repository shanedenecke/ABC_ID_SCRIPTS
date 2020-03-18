#!/usr/bin/env bash
H=~/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H

#### CAFE
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees

source ./ABC_ID_SCRIPTS/ABC_Ultrametric_tree_generate.sh