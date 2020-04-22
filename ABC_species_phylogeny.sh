cd $H 
cd CAFE

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Hemiptera_species_taxid.txt -ortho_algo Orthofinder -outgroups "DroMel BomMor ApiMel AedAeg " -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Lepidoptera_species_taxid.txt -ortho_algo Orthofinder -outgroups "MyzPer" -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Diptera_species_taxid.txt -ortho_algo OrthoDB -outgroups "BomMor" -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Arthropod_species_taxid.txt -ortho_algo Orthofinder -outgroups "CaeEle" -threads $THREADS

cp ./*/rax_output/RAxML_bipartitions.*.nwk ./clean_raxml_trees/