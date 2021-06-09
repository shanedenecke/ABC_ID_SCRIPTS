#!/usr/bin/env bash

### 1) Set environemntal variables for pipeline
H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=80
THREADS=16
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
 ~/Applications/ABC_scan/ABC_scan.sh -proteome $i -outdir ./ABC_search -threads $THREADS -minlen 250
done

###### 5)Make figures and tables from counts
Rscript ./ABC_ID_SCRIPTS/ABC_post_process.R

###### 6) CAFE
mkdir CAFE ./CAFE/clean_raxml_trees ./CAFE/taxid_lists ./CAFE/outputs

### Prepare CAFE
./ABC_ID_SCRIPTS/Species_phylogeny_prep.py
#for i in ./CAFE/taxid_lists/*.txt;do ~/Applications/Custom_Applications/Species_phylogeny.sh -taxid $i -threads $THREADS -outdir ./CAFE -maxseqs 1000; done
cp ./GENERAL_REFERENCE/CAFE/ultrametric_tree_backup/*.support ./CAFE/clean_raxml_trees/
Rscript ./ABC_ID_SCRIPTS/ABC_CAFE5_prep.R  

### Run CAFE
for i in ./CAFE/CAFE_tables/*.tsv; do  
b=$(echo $(basename $i) | sed 's/_ABC_CAFE_table.tsv//g')
cafexp -i $i -o ./CAFE/outputs/$b -t ./CAFE/clean_raxml_trees/$b'_tree_ultrametric.nwk' 
done 

##Make CAFE figures
Rscript ./ABC_ID_SCRIPTS/ABC_CAFE5_figures.R
  

  

### 7) Make phylogeny for relevant species 
./ABC_ID_SCRIPTS/ABC_phylo.sh


