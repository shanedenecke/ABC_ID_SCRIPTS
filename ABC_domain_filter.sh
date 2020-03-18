#!/usr/bin/env bash
cd $H

mkdir Filter

### get all fasta files
#z=$(($THREADS / 4))
cat ./preliminary_ABC/proteomes/* > ./Filter/ABC_preliminary_total.faa
cat ./Db_build_temp/Only_ABCs.faa* >> ./Filter/ABC_preliminary_total.faa

### Run IP scan and extract domains 
#~/interproscan-5.30-69.0/interproscan.sh -cpu $z -appl pfam -i ./Filter/ABC_preliminary_total.faa -o ./Filter/IPSCAN.tsv -f TSV
hmmsearch --domtblout ./Filter/HMM_PF00005_output.tsv ./ABC_REF/Model_ABC_sets/ABC_tran.hmm ./Filter/ABC_preliminary_total.faa > ./Filter/hmm_junk.txt
cat ./Filter/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > ./Filter/HMM_PF00005_output_clean.tsv

#~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./Filter/ABC_preliminary_total.faa -table ./Filter/HMM_PF00005_output.tsv > ./Filter/Only_ABCs_NBD.faa

#### Run R script to do preliminary sort 
Rscript ./ABC_ID_SCRIPTS/ABC_Domain_Filter_Prelim.R $QUAL_THRESH

### Make unrooted tree for all unsorted proteins 
#mkdir ./Filter/Unsorted_clean
#x=./Filter/Unsorted_clean/ABC_unsorted_total
#grep -A 1 '_Unsorted' ./Filter/ABC_preliminary_total.faa | sed '/--/d' > $x.faa
#grep -E -A 1 "DroMel" ./Db_build_temp/Only_ABCs.faa | sed '/--/d' >> $x.faa
#~/Applications/Custom_Applications/fasta_seq_len_filter.py $x.faa 300 10000 > $x'_length.faa'
#~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta $x'_length.faa' -table ./Filter/HMM_PF00005_output.tsv > $x'_length_NBD.faa'
#grep -A 1 -f ./Filter/preliminary/Quality_threshold_species.txt $x'_length_NBD.faa' | sed '/--/d' > $x'_length_NBD_species.faa'
#cat $x'_length_NBD_species.faa' | sed 's/\:/_/g' | sed 's/\-/_/g' | sed 's/\./_/g' > $x'_length_NBD_quality_clean.faa'
#~/Applications/Custom_Applications/fasta_remove.py $x'_length_NBD_quality_clean.faa' ./ABC_REF/non_model_proteomes/bad_ABCs.txt > $x'_length_NBD_quality_clean_good.faa'

#y=$x'_length_NBD_quality_clean_good.faa'
#mafft --thread $THREADS $y > $y.aln
#~/Applications/trimAl/source/trimal -automated1 -in $y'.aln' -out $y'.aln.trimm' 
#~/Applications/Custom_Applications/fasta_2_phylip.sh $y'.aln.trimm' > $y'.aln.trimm.phy'
#Rscript ~/Applications/Custom_Applications/Phylip_duplicate.R $y'.aln.trimm.phy' > $y'_clean.phy'

### run tree
#nohup ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 200 -T $THREADS -m PROTGAMMAAUTO -s $y'_clean.phy' -n 'Unsorted.nwk' -w $H/Filter/Unsorted_clean &
# cp 

Rscript ./ABC_ID_SCRIPTS/ABC_Domain_Filter_include_unsort.R