#!/usr/bin/env bash

### Make proteome output directory
mkdir proteome_clean
mkdir ./proteome_clean/clean_fasta

### Get list of taxid codes to be used for your species. This uses the $SPEC variable which is $H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv file. Subset this file if you only want to do certain species
### Get rid of model proteomes as they will be added later
cat $SPEC | awk '{print $2 "\t" $1}'| tail -n +2 | grep -v "DroMel" | grep -v "HomSap" | grep -v "TriCas" | grep -v "TetUrt" > ./proteome_clean/target_species.tsv

### Perform loop to go through taxid list and get those OrhtoDB files. Rename the files with orthoDB codes
while IFS=$'\t' read -r taxid abbrev
do
  if [ -e ./GENERAL_REFERENCE/non_model_proteomes/orthoDB_fasta/$taxid'_0.fs' ]; then
    #echo $abbrev
    cp ./GENERAL_REFERENCE/non_model_proteomes/orthoDB_fasta/$taxid'_0'* ./proteome_clean/clean_fasta/"$abbrev".fasta
    grep -e "^$taxid" ./GENERAL_REFERENCE/non_model_proteomes/keys/odb10v0_genes.tab | cut -f 1,3 > ./proteome_clean/clean_fasta/"$abbrev"_OrthoDB_key.txt
    echo -e "code,name" | cat - ./proteome_clean/clean_fasta/"$abbrev"_OrthoDB_key.txt | perl -pe 's/\t/,/g' > ./proteome_clean/clean_fasta/"$abbrev"_OrthoDB_key_DICT.txt
    ~/Applications/Custom_Applications/fasta_rename_fast_only_exact.py ./proteome_clean/clean_fasta/"$abbrev".fasta ./proteome_clean/clean_fasta/"$abbrev"_OrthoDB_key_DICT.txt > ./proteome_clean/clean_fasta/"$abbrev"_unigene.faa
  else 
    echo $abbrev' not found in OrthoDB'
  fi
done < ./proteome_clean/target_species.tsv


#### Remove empty or duplicate sequences
#find ./proteome_clean/clean_fasta/ -size  0 -print0 | xargs -0 rm --
for i in ./proteome_clean/clean_fasta/*_unigene.faa
do
sed 's/\.//g' $i | sed 's/\*//g' | sed 's/ /_/g' | sed 's/\\//g' | sed 's/\///g' | tr '[:lower:]' '[:upper:]' > temp.fasta
awk '/^>/{id=$0;getline;arr[id]=$0}END{for(id in arr)printf("%s\n%s\n",id,arr[id])}' temp.fasta > $i
done
rm temp.fasta

#### makeblastDB for each proteome
for i in ./proteome_clean/clean_fasta/*_unigene.faa
do
makeblastdb -in $i -parse_seqids -dbtype prot
done

### report number of good and total proteomes
echo 'Total numberof proteomes is ' $(ls ./proteome_clean/clean_fasta/* | grep -E "*unigene.faa$" | wc -l)
echo 'Numberof GOOD proteomes is ' $(ls ./proteome_clean/clean_fasta/* | grep -E "*psq$" | wc -l)


mv ./proteome_clean/clean_fasta/*_unigene.faa ./proteomes/
cat ./proteome_clean/target_species.tsv | cut -f 2 | while read i;do
cp ./GENERAL_REFERENCE/non_model_proteomes/non_orthoDB_fasta/$i*.faa ./proteomes/
done

rm -rf ./proteome_clean/

for i in ./proteomes/*; do fasta_clean.sh -proteome $i; done
