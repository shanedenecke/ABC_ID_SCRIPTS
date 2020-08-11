#!/usr/bin/env bash
H=/mnt/disk/shane/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.3
THREADS=10

cd $H 

cd CAFE

~/Applications/Custom_Applications/Species_phylogeny.sh -threads 10 -taxid_codes ../ABC_REF/CAFE/Lepidoptera_species.txt -ortho_algo Orthofinder -outgroups "ApiMel" -threads $THREADS