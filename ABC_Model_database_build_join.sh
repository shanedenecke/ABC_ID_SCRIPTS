#!/usr/bin/env bash

cd $H
mkdir Db_build_temp

#### Concatonate all model proteomes and dictionaries
cat ./ABC_REF/Model_ABC_sets/Raw_proteomes/*.faa > ./Db_build_temp/combined_proteom.faa
cat ./ABC_REF/Model_ABC_sets/Dicts/*.csv | grep -v 'code' > ./Db_build_temp/combined_dictionary_nohead.csv
echo -e 'code,name' | cat - ./Db_build_temp/combined_dictionary_nohead.csv > ./Db_build_temp/combined_dictionary.csv

#### Rename master combined fasta file: Marks for ABC transporters 
$H/ABC_ID_SCRIPTS/general_scripts/fasta_rename.py ./Db_build_temp/combined_proteom.faa ./Db_build_temp/combined_dictionary.csv > ./Db_build_temp/combined_marked_proteome.faa
makeblastdb -in ./Db_build_temp/combined_marked_proteome.faa -parse_seqids -dbtype prot


###### Extract all model ABC sequences into separate family fasta files, align them and then create HMM database
grep -A 1 -E "__ABC[A-Z]+" ./Db_build_temp/combined_marked_proteome.faa | sed '/--/d' > ./Db_build_temp/Only_ABCs.faa
/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i Db_build_temp/Only_ABCs.faa -o Db_build_temp/IPSCAN -f TSV
/data2/shane/Applications/custom/ip_domain_extract.py ./Db_build_temp/Only_ABCs.faa ./Db_build_temp/IPSCAN "PF00005" > ./Db_build_temp/Only_ABCs_NBD.faa

mafft --thread $THREADS ./Db_build_temp/Only_ABCs.faa > ./Db_build_temp/Only_ABCs.aln
$H'/ABC_ID_SCRIPTS/general_scripts/trimAl/source/trimal' -in ./Db_build_temp/Only_ABCs.aln -out ./Db_build_temp/Only_ABCs.trimmed
hmmbuild ./Db_build_temp/Only_ABCs.hmm ./Db_build_temp/Only_ABCs.trimmed


mkdir model_database
mkdir ./model_database/HMM_databases
mkdir ./model_database/marked_proteome
mkdir ./model_database/model_NBDs

mv ./Db_build_temp/*.hmm ./model_database/HMM_databases/
mv ./Db_build_temp/combined_marked_proteome.fa* ./model_database/marked_proteome/
mv ./Db_build_temp/*NBD* ./model_database/model_NBDs/


#cat ./ABC_REF/Input_files/ABC_families.txt | while read i
#do
#  b='./Db_build_temp/Model_'$(echo $i | sed 's/_//g')'_NBD.faa'
#  grep -A 1 $i ./Db_build_temp/Only_ABCs_NBD.faa | sed '/--/d' > $b
#  mafft --thread $THREADS $b > $b.aln
#  $H'/ABC_ID_SCRIPTS/general_scripts/trimAl/source/trimal' -in $b.aln -out $b.trimmed
#  hmmbuild $b.hmm $b.trimmed
#done

