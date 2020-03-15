#!/usr/bin/env bash
H='~/Transporter_ID/ABC_ID'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H

###############1) proteome prepare
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh

