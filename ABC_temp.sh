#!/usr/bin/env bash
H=/mnt/disk/shane/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.3
THREADS=14

cd $H 

source ./ABC_ID_SCRIPTS/ABC_phylo.sh