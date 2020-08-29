#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=80
THREADS=4
cd $H 


### 3) Prepare proteomes
mkdir proteomes
source ./ABC_ID_SCRIPTS/ABC_proteome_prepare.sh
  
###### 2) Build database materials
source ./ABC_ID_SCRIPTS/ABC_Model_database_build.sh


##### 3) Run BUSCO for quality assessment
source ./ABC_ID_SCRIPTS/ABC_BUSCO.sh

###### 4) Search proteomes
for i in ./proteomes/*.faa; do
 ./ABC_ID_SCRIPTS/ABC_search/ABC_search.sh -proteome $i -outdir ./ABC_search -threads $THREADS -minlen 250
done

###### 5)Make figures and tables from counts
Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R

###### 6) CAFE
mkdir CAFE
mkdir ./CAFE/clean_raxml_trees
mkdir ./CAFE/taxid_lists

./ABC_ID_SCRIPTS/Species_phylogeny_prep.py
#for i in ./CAFE/taxid_lists/*.txt;do ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid_codes $i -threads $THREADS -outdir ./CAFE -maxseqs 1000; done
#nohup ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid ./CAFE/taxid_lists/Hemimetabola_taxid_codes.txt -threads $THREADS -outdir ./CAFE -maxseqs 500 &
#nohup ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid ./CAFE/taxid_lists/Coleoptera_taxid_codes.txt -threads $THREADS -outdir ./CAFE -maxseqs 500 &
#nohup ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid ./taxid_lists/Lepidoptera_taxid_codes.txt -threads $THREADS -maxseqs 500 &
#nohup ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid ./CAFE/taxid_lists/Arachnid_taxid_codes.txt -threads $THREADS -outdir ./CAFE -maxseqs 500 &
cp ./GENERAL_REFERENCE/CAFE/ultrametric_tree_backup/*.support ./CAFE/clean_raxml_trees/

#Rscript ./ABC_ID_SCRIPTS/ABC_CAFE_prep.R  
#source ./ABC_ID_SCRIPTS/ABC_CAFE5_run_full.sh
#Rscript ./ABC_ID_SCRIPTS/ABC_CAFE5_figures.R
  

  

### 7) Make phylogeny for relevant species 

#source ./ABC_ID_SCRIPTS/ABC_phylo.sh
#Rscript ~/Applications/Custom_Applications/ggtree_clean_phylogeny.R #### NEED TO GET ALL TREES IN ONE FOLDER
#cp ./GENERAL_REFERENCE/phylo_premade/* ./phylo/clean_trees/  


