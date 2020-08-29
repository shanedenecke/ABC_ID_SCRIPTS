#!/usr/bin/env bash


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"




### add help 
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search SLC transporters in non-model arthropods
  
  Arguments:
  
  -target: Path to folder containing one or more Arthropod proteomes that you wish to search. Protoemes should be labled witht a 6 letter abbreviation followed by '_unigene.faa' e.g. DroMel_unigene.faa for Drosophila melanogaster 
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: They will be put into relevant subdiretories automatically
  -hmm_profile: A HMM profile to search with
  -database: protein database for reciprocal blast
  -print: Do you want to print the final output to the terminal
  CRIPTS/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
  "
  exit 0
fi

### add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
    -hmm_profile) HMM_PROFILE="$2"; shift 2;;
    -database) DATABASE="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -print) PRINT="$2"; shift 2;;
    -motif) MOTIF="$2"; shift 2;;
    -minlen) MINLEN="$2"; shift 2;;
    -motif) MOTIF="$2"; shift 2;;
    -domain_filter) DOMAIN_FILTER="$2"; shift 2;;
  esac  
done


### Debugging
#cd /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
#PROTEOME=./proteomes/DanPle_unigene.faa
#SCRIPT_DIR='/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/ABC_ID_SCRIPTS/ABC_search'


if [ -z $HMM_PROFILE ];then HMM_PROFILE=$SCRIPT_DIR/ABC_tran.hmm; fi
if [ -z $MOTIF ];then MOTIF=$SCRIPT_DIR/ABC_meme_clean.txt; fi
if [ -z $DATABASE ];then DATABASE=$SCRIPT_DIR/combined_marked_proteome.faa; fi
if [ -z $OUTDIR ];then OUTDIR="./ABC_search"; fi
if [ -z $PRINT ];then PRINT="NO"; fi
if [ -z $THREADS ];then THREADS=4; fi
if [ -z $MINLEN ];then MINLEN=250; fi
if [ -z $MOTIF ];then MOTIF=NA; fi
if [ -z $DOMAIN_FILTER ];then DOMAIN_FILTER=NA; fi

#echo $PROTEOME
#echo $THREADS
#echo $HMM_PROFILE
#echo $DATABASE


### establish basename for proteome
bname=$(basename $PROTEOME | sed 's/_unigene.faa//g')
mkdir -p $OUTDIR
mkdir $OUTDIR/$bname
echo ' Now searching the '$bname' proteome'


### perform HMM search
hmmsearch --cpu $THREADS --notextw -E 100 $HMM_PROFILE $PROTEOME | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > $OUTDIR/$bname/total_ABC_HMM.table


### clean HMM search and output fasta
cut -f 9 $OUTDIR/$bname/total_ABC_HMM.table | sed 's/\s+//g' | ~/Applications/Custom_Applications/unigene_fa_sub.sh $PROTEOME - > $OUTDIR/$bname/total_ABC_pre_blast.faa

echo 'BLAST AWAY'
## run reciprocal blast
blastp -query $OUTDIR/$bname/total_ABC_pre_blast.faa -db $DATABASE -outfmt "6 qseqid sseqid pident evalue qlen sstart send" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 -num_threads $THREADS > $OUTDIR/$bname/total_ABC_recip_blast.tsv

### sort into families 
Rscript $SCRIPT_DIR/ABC_Family_Sort.R $OUTDIR/$bname/total_ABC_recip_blast.tsv > $OUTDIR/$bname/total_ABC_dict.csv

## get preliminary fasta files 
~/Applications/Custom_Applications/fasta_rename.py $PROTEOME $OUTDIR/$bname/total_ABC_dict.csv > $OUTDIR/$bname/total_ABC_marked_proteome.faa
grep -A 1 "__ABC" $OUTDIR/$bname/total_ABC_marked_proteome.faa | sed '/--/d' > $OUTDIR/$bname/Preliminary_ABCs.faa

rm $OUTDIR/$bname/total_ABC_marked_proteome.faa

## Get domain content
hmmsearch --domtblout $OUTDIR/$bname/HMM_PF00005_output.tsv $HMM_PROFILE $OUTDIR/$bname/Preliminary_ABCs.faa > $OUTDIR/$bname/hmm_junk.txt
cat $OUTDIR/$bname/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > $OUTDIR/$bname/HMM_PF00005_output_clean.tsv

## Get motif content
mast -hit_list $MOTIF $OUTDIR/$bname/Preliminary_ABCs.faa -notext -nohtml -nostatus | tail -n +3 | head -n -1 | sed -E 's/\s+/___/g' | sed 's/___/\t/g' | cut -f 1,3 > $OUTDIR/$bname/Mem_motif.tsv


### Make summary table
#$SCRIPT_DIR/ABC_domain_filter.R $(readlink -f $OUTDIR/$bname)
$SCRIPT_DIR/ABC_table_summarize.R --minlen 250 --motif $MOTIF --domain_filter $DOMAIN_FILTER --indir $(readlink -f $OUTDIR/$bname)


mkdir $OUTDIR/$bname/junk
mv $OUTDIR/$bname/total_ABC_HMM.table $OUTDIR/$bname/total_ABC_pre_blast.faa $OUTDIR/$bname/ABC_filtered_out.csv $OUTDIR/$bname/HMM_PF00005_output.tsv $OUTDIR/$bname/HMM_PF00005_output_clean.tsv $OUTDIR/$bname/Mem_motif.tsv $OUTDIR/$bname/junk
mv $OUTDIR/$bname/hmm_junk.txt $OUTDIR/$bname/Preliminary_ABCs.faa $OUTDIR/$bname/total_ABC_dict.csv $OUTDIR/$bname/total_ABC_recip_blast.tsv $OUTDIR/$bname/junk

if [ $PRINT == "YES" ]; then cat $OUTDIR/$bname/Final_ABC_table.tsv; fi


