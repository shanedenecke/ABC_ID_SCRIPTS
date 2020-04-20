cd $H 
cd CAFE

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Hemi_new.txt -ortho_algo Orthofinder -outgroups "DroMel BomMor ApiMel AedAeg " -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Lep_new.txt -ortho_algo Orthofinder -outgroups "MyzPer" -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Dipt_new.txt -ortho_algo OrthoDB -outgroups "None" -threads $THREADS

~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes ../ABC_REF/CAFE/Arth_new3.txt -ortho_algo Orthofinder -outgroups "CaeEle" -threads $THREADS

