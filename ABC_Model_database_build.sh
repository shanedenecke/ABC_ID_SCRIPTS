#!/usr/bin/env bash

cd $H
mkdir Db_build_temp

#### Concatonate all model proteomes and dictionaries
cat ./GENERAL_REFERENCE/Model_ABC_sets/Raw_proteomes/*.faa > ./Db_build_temp/combined_proteom.faa
cat ./GENERAL_REFERENCE/Model_ABC_sets/Dicts/*.csv | grep -v 'code' > ./Db_build_temp/combined_dictionary_nohead.csv
echo -e 'code,name' | cat - ./Db_build_temp/combined_dictionary_nohead.csv > ./Db_build_temp/combined_dictionary.csv

#### Rename master combined fasta file: Marks for ABC transporters 
~/Applications/Custom_Applications/fasta_rename.py ./Db_build_temp/combined_proteom.faa ./Db_build_temp/combined_dictionary.csv > ./Db_build_temp/combined_marked_proteome.faa
makeblastdb -in ./Db_build_temp/combined_marked_proteome.faa -parse_seqids -dbtype prot


###### Extract all model ABC sequences, align them and then create HMM database
grep -A 1 -E "__ABC[A-Z]+" ./Db_build_temp/combined_marked_proteome.faa | sed '/--/d' > ./Db_build_temp/Only_ABCs.faa
#~/Applications/interproscan-5.30-69.0/interproscan.sh -appl pfam -i Db_build_temp/Only_ABCs.faa -o Db_build_temp/IPSCAN -f TSV
#~/Applications/Custom_Applications/ip_domain_extract.py ./Db_build_temp/Only_ABCs.faa ./Db_build_temp/IPSCAN "PF00005" > ./Db_build_temp/Only_ABCs_NBD.faa

hmmsearch --domtblout ./Db_build_temp/PF0005_HMM_output.tsv ./GENERAL_REFERENCE/Model_ABC_sets/ABC_tran.hmm ./Db_build_temp/Only_ABCs.faa > ./Db_build_temp/hmm_junk.txt
~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./Db_build_temp/Only_ABCs.faa -table ./Db_build_temp/PF0005_HMM_output.tsv > ./Db_build_temp/Only_ABCs_NBD.faa

mafft-linsi --thread $THREADS ./Db_build_temp/Only_ABCs.faa > ./Db_build_temp/Only_ABCs.aln
~/Applications/trimal/source/trimal -automated1 -in ./Db_build_temp/Only_ABCs.aln -out ./Db_build_temp/Only_ABCs.trimmed
hmmbuild ./Db_build_temp/Only_ABCs.hmm ./Db_build_temp/Only_ABCs.trimmed


mkdir model_database
mkdir ./model_database/HMM_databases
mkdir ./model_database/marked_proteome
mkdir ./model_database/model_NBDs

mv ./Db_build_temp/*.hmm ./model_database/HMM_databases/
mv ./Db_build_temp/combined_marked_proteome.fa* ./model_database/marked_proteome/
mv ./Db_build_temp/*NBD* ./model_database/model_NBDs/
mv ./Db_build_temp/Only_ABCs.faa ./model_database/Model_organism_ABCs.faa
mv ./Db_build_temp ./intermediate
