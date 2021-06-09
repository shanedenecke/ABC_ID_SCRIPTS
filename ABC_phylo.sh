#!/usr/bin/env bash

H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
THREADS=8
cd $H 



mkdir phylo
mkdir ./phylo/clean_trees

cd phylo
########## Make ABCH tree specifically
#grep -E -A 1 'NezVir|DakVit|MyzPer|BemTab|NilLug|DapPul|TetUrt|RhoPro' ../Final_outputs/combined_files/All_ABCs.faa | grep -A 1 'ABCH' | sed '/--/d' > ./ABCH_total.faa
#remove_duplicates.sh ABCH_total.faa > ABCH_total_unique.faa
#nohup align_and_tree.sh -proteins /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/phylo/ABCH_total_unique.faa -threads $THREADS & 
#cp ./ABCH/*support 

########## Make ABCBF tree specifically
grep -E -A 1 'HelArm|BomMor|AmyTra|PapPol|PieRap|DanPle|NezVir' ../Final_outputs/combined_files/All_ABCs.faa | grep -A 1 'ABCBF' | sed '/--/d' > ./ABCBF_total.faa
align_and_tree.sh -proteins /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/phylo/ABCBF_total.faa -threads $THREADS
cp ./ABCBF/*support ./clean_trees

########## Make ABCCC tree specifically
grep -E -A 1 'LepDec|DenPon|TriCas|OntTau|DroMel' ../Final_outputs/combined_files/All_ABCs.faa | grep -A 1 'ABCC' | sed '/--/d' > ./ABCC_total.faa
align_and_tree.sh -proteins /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/phylo/ABCC_total.faa -threads $THREADS
cp ./ABCC/*support ./clean_trees










