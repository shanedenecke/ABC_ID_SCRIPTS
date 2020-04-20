#!/usr/bin/env bash
H=/mnt/disk/shane/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.3
THREADS=10

cd $H 
cd CAFE
#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Hemi_new.txt -ortho_algo Orthofinder -outgroups "DroMel BomMor ApiMel AedAeg " -threads $THREADS
#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Lep_new.txt -ortho_algo Orthofinder -outgroups "MyzPer" -threads $THREADS
#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Dipt_new.txt -ortho_algo OrthoDB -outgroups "None" -threads $THREADS
~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Arth_new3.txt -ortho_algo Orthofinder -outgroups "CaeEle" -threads $THREADS