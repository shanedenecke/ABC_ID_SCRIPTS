#!/usr/bin/env bash
#$1 == 'proteome that you want to search'
#i=./proteomes/ApiMel_unigene.faa

cd $H
bname=$(basename $1 | sed 's/_unigene.faa//g')
echo $bname
#bname=$(basename $1 | sed 's/_unigene.faa//g')

mkdir ABC_search/$bname

### perform HMM search
#hmmsearch --cpu $THREADS --notextw -E 20 ./model_database/HMM_databases/Only_ABCs.hmm $1 | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > ./ABC_search/$bname/total_ABC_HMM.table
hmmsearch --cpu $THREADS --notextw -E 20 ./ABC_REF/Model_ABC_sets/ABC_tran.hmm $1 | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > ./ABC_search/$bname/total_ABC_HMM.table
#hmmsearch --notextw -E 20 ./model_database/HMM_databases/Only_ABCs.hmm $i | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > ./ABC_search/$bname/total_ABC_HMM.table

### clean HMM search
cut -f 9 ./ABC_search/$bname/total_ABC_HMM.table | sed 's/\s+//g' \
| ~/Applications/Custom_Applications/unigene_fa_sub.sh $1 - > ./ABC_search/$bname/total_ABC_pre_blast.faa

#cut -f 9 ./ABC_search/$bname/total_ABC_HMM.table | sed 's/\s+//g' \
#| ~/Applications/Custom_Applications/unigene_fa_sub.sh $i - > ./ABC_search/$bname/total_ABC_pre_blast.faa

echo 'BLAST AWAY'
## run reciprocal blast
blastp -query ./ABC_search/$bname/total_ABC_pre_blast.faa -db ./model_database/marked_proteome/combined_marked_proteome.faa -outfmt "6 qseqid sseqid pident evalue qlen sstart send" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 -num_threads $THREADS > ./ABC_search/$bname/total_ABC_recip_blast.tsv

### sort into families 
Rscript ./ABC_ID_SCRIPTS/ABC_Family_Sort.R ./ABC_search/$bname/total_ABC_recip_blast.tsv > ./ABC_search/$bname/total_ABC_dict.csv

## get preliminary fasta files 
~/Applications/Custom_Applications/fasta_rename.py ./proteomes/$bname'_unigene.faa' ./ABC_search/$bname/total_ABC_dict.csv > ./ABC_search/$bname/total_ABC_marked_proteome.faa
grep -A 1 "__ABC" ./ABC_search/$bname/total_ABC_marked_proteome.faa | sed '/--/d' > ./ABC_search/$bname/Preliminary_ABCs.faa
rm ./ABC_search/$bname/total_ABC_marked_proteome.faa


mv ./ABC_search/$bname/total_ABC_dict.csv ./preliminary_ABC/dicts/$bname'_preliminary_ABC_dict.csv'
mv ./ABC_search/$bname/Preliminary_ABCs.faa ./preliminary_ABC/proteomes/$bname'_preliminary_ABC.faa'
 