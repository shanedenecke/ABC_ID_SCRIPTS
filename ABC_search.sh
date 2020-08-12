#!/usr/bin/env bash



### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -proteomes: Path to folder containing one or more Arthropod proteomes that you wish to search. Protoemes should be labled witht a 6 letter abbreviation followed by '_unigene.faa' e.g. DroMel_unigene.faa for Drosophila melanogaster 
  -busco_thresh: Threshold that you wish to set for BUSCO completeness scores. e.g. 75 means that only proteomes with >75 BUSCO score will be considered for the analysis
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
  -metadata: A tab separated table containign the 'Species_name' and the 'abbreviation' of the targeted species
  example
  ./SLC_ID_SCRIPTS/SLC_id_standalone/SLC_id_scripts/SLC_id.sh -proteomes $H/proteomes -busco_thresh 75 -threads $THREADS -outdir $H -metadata /mnt/disk/shane/Transporter_ID/SLC_id_pipeline/SLC_ID_SCRIPTS/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -target) PROTEOME="$2"; shift 2;;
    -hmm_profile) HMM_PROFILE="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
  esac
done


### Debugging
#PROTEOME=./proteomes/HelVir_unigene.faa
#HMM_PROFILE=./model_database/HMM_databases/Only_ABCs.hmm
#THREADS=14
#OUTDIR=./ABC_search




### establish basename for proteome
bname=$(basename $PROTEOME | sed 's/_unigene.faa//g')
mkdir $OUTDIR/$bname
echo ' Now searching the '$bname' proteome'


### perform HMM search
hmmsearch --cpu $THREADS --notextw -E 20 $HMM_PROFILE $PROTEOME | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > ./ABC_search/$bname/total_ABC_HMM.table


### clean HMM search and output fasta
cut -f 9 ./ABC_search/$bname/total_ABC_HMM.table | sed 's/\s+//g' | ~/Applications/Custom_Applications/unigene_fa_sub.sh $PROTEOME - > ./ABC_search/$bname/total_ABC_pre_blast.faa

echo 'BLAST AWAY'
## run reciprocal blast
blastp -query ./ABC_search/$bname/total_ABC_pre_blast.faa -db ./model_database/marked_proteome/combined_marked_proteome.faa -outfmt "6 qseqid sseqid pident evalue qlen sstart send" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 -num_threads $THREADS > ./ABC_search/$bname/total_ABC_recip_blast.tsv

### sort into families 
Rscript ./ABC_ID_SCRIPTS/ABC_Family_Sort.R ./ABC_search/$bname/total_ABC_recip_blast.tsv > ./ABC_search/$bname/total_ABC_dict.csv

## get preliminary fasta files 
~/Applications/Custom_Applications/fasta_rename.py ./proteomes/$bname'_unigene.faa' ./ABC_search/$bname/total_ABC_dict.csv > ./ABC_search/$bname/total_ABC_marked_proteome.faa
grep -A 1 "__ABC" ./ABC_search/$bname/total_ABC_marked_proteome.faa | sed '/--/d' > ./ABC_search/$bname/Preliminary_ABCs.faa
rm ./ABC_search/$bname/total_ABC_marked_proteome.faa

## Filter based on domain content
hmmsearch --domtblout ./ABC_search/$bname/HMM_PF00005_output.tsv ./GENERAL_REFERENCE/Model_ABC_sets/ABC_tran.hmm ./ABC_search/$bname/Preliminary_ABCs.faa > ./ABC_search/$bname/hmm_junk.txt
cat ./ABC_search/$bname/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > ./ABC_search/$bname/HMM_PF00005_output_clean.tsv


./ABC_ID_SCRIPTS/ABC_domain_filter.R $(readlink -f ./ABC_search/$bname)
