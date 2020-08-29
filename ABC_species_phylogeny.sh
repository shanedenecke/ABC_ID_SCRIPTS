#!/usr/bin/env bash

cd /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
cd CAFE

for i in ./taxid_lists/*.txt;do ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes $i -ortho_algo Orthofinder -outgroups "None" -threads 10; done


~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ./taxid_lists/Hemimetabola_taxid_codes.txt -ortho_algo Orthofinder -outgroups "None" -threads 10; done


#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/CAFE/Hemiptera_species.txt -ortho_algo Orthofinder -outgroups "DroMel BomMor ApiMel AedAeg " -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/CAFE/Lepidoptera_species.txt -ortho_algo Orthofinder -outgroups "ApiMel" -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../GENERAL_REFERENCE/CAFE/Diptera_species.txt -ortho_algo OrthoDB -outgroups "BomMor" -threads $THREADS

#~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ./GENERAL_REFERENCE/CAFE/ArachInsect_species.txt -ortho_algo Orthofinder -outgroups "None" -threads 10

#cp ./*/rax_output/RAxML_bipartitions.*.nwk ./clean_raxml_trees/
