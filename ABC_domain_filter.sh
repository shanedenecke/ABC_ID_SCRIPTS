#!/usr/bin/env bash
cd $H





mkdir Filter

### get all fasta files
cat ./preliminary_ABC/proteomes/* > ./Filter/ABC_preliminary_total.faa
cat ./Db_build_temp/Only_ABCs.faa* >> ./Filter/ABC_preliminary_total.faa

### Run IP scan and extract domains 
/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i ./Filter/ABC_preliminary_total.faa -o ./Filter/IPSCAN.tsv -f TSV

Rscript ./ABC_ID_SCRIPTS/ABC_Domain_Filter.R $QUAL_THRESH


### Make unrooted tree for all unsorted proteins 
#grep -A 1 '_Unsorted' ./Filter/ABC_preliminary_total.faa | sed '/--/d' > ./Filter/ABC_unsorted_total.faa
#/data2/shane/Applications/custom/ip_domain_extract.py ./Filter/ABC_unsorted_total.faa ./Filter/IPSCAN.tsv "PF00005" | sed 's/:/_/g' > ./Filter/ABC_unsorted_total_NBD.faa 
#cat ./model_database/model_NBDs/Only_ABCs_NBD.faa | grep -E -A 1 'DroMel' >> ./Filter/ABC_unsorted_total_NBD.faa 
#grep -f ./Filter/Quality_cutoff_species.txt ./Filter/ABC_unsorted_total_NBD.faa > ./Filter/ABC_unsorted_total_NBD_filter.faa
#cd Filter
#nohup /data2/shane/Applications/custom/align_and_tree.sh ./ABC_unsorted_total_NBD.faa $THREADS &
#cd ..

