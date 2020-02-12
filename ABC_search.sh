#!/usr/bin/env bash
#$1 == 'proteome that you want to search'
#A=./proteomes/MyzPer_unigene.faa

cd $H
bname=$(basename $1 | sed 's/_unigene.faa//g')
#bname=$(basename $A | sed 's/_unigene.faa//g')

mkdir ABC_search
mkdir ABC_search/$bname
#cp ./proteomes/$bname'_unigene.faa' ./ABC_search/$bname/

for i in ./model_database/HMM_databases/*; do
  database=$(echo $(basename $i | perl -pe 's/Model_(ABC.+).faa.hmm/$1/g'))
  famname=$(echo $database | sed 's/_NBD//g')
  #hmmsearch --notextw -E 20 $i $1 > ./ABC_search/$bname/$database.hmmoutput
  
  hmmsearch --notextw -E 20 $i $1 | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' \
  | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" \
  | cut -f 2- > ./ABC_search/$bname/$database.table
  
  cut -f 9 ./ABC_search/$bname/$database.table | sed 's/\s+//g' \
  | /data2/shane/Applications/custom/unigene_fa_sub.sh $1 - > ./ABC_search/$bname/$database.faa
  
  blastp -query ./ABC_search/$bname/$database.faa -db ./model_database/marked_proteome/combined_marked_proteome.faa\
  -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 -num_threads $THREADS > ./ABC_search/$bname/$database'_recip_blast.tsv'
  
  Rscript $H/ABC_ID_SCRIPTS/ABC_Family_Sort.R ./ABC_search/$bname/$database'_recip_blast.tsv' > ./ABC_search/$bname/$database'_ABC_dict.csv'
  
  $H/ABC_ID_SCRIPTS/general_scripts/fasta_rename.py ./proteomes/$bname'_unigene.faa' ./ABC_search/$bname/$database'_ABC_dict.csv' > ./ABC_search/$bname/$database'_ABC_marked_proteome.faa'
  
  grep -A 1 "__ABC" ./ABC_search/$bname/$database'_ABC_marked_proteome.faa' | sed '/--/d' > ./ABC_search/$bname/$database'_preliminary_ABCs.faa'
  
  rm ./ABC_search/$bname/$database'_ABC_marked_proteome.faa'
done


cat ./ABC_search/$bname/*'_preliminary_ABCs.faa'  | sed '/name/d' > ./ABC_search/$bname/$bname'_total_preliminary_ABCs.faa'
mv ./ABC_search/$bname/$bname'_total_preliminary_ABCs.faa' preliminary_ABC/proteomes

cat ./ABC_search/$bname/*'_ABC_dict.csv' | sed '/name/d'  | sed 1i"code,name" > ./ABC_search/$bname/$bname'_total_preliminary_ABCs_dict.csv'
mv ./ABC_search/$bname/$bname'_total_preliminary_ABCs.faa' preliminary_ABC/dicts