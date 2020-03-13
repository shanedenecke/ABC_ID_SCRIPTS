#!/usr/bin/env bash
cd $H

mkdir Filter

### get all fasta files
z=$(($THREADS / 4))
cat ./preliminary_ABC/proteomes/* > ./Filter/ABC_preliminary_total.faa
cat ./Db_build_temp/Only_ABCs.faa* >> ./Filter/ABC_preliminary_total.faa

### Run IP scan and extract domains 
/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -cpu $z -appl pfam -i ./Filter/ABC_preliminary_total.faa -o ./Filter/IPSCAN.tsv -f TSV

#### Run R script to do preliminary sort 
Rscript ./ABC_ID_SCRIPTS/ABC_Domain_Filter_Prelim.R $QUAL_THRESH

### Make unrooted tree for all unsorted proteins 
mkdir ./Filter/Unsorted_clean
x=./Filter/Unsorted_clean/ABC_unsorted_total
grep -A 1 '_Unsorted' ./Filter/ABC_preliminary_total.faa | sed '/--/d' > $x.faa
grep -E -A 1 "DroMel|HomSap" ./Db_build_temp/Only_ABCs.faa | sed '/--/d' >> $x.faa
/data2/shane/Applications/custom/fasta_seq_len_filter.py $x.faa 100 10000 > $x'_length.faa'
/data2/shane/Applications/custom/ip_domain_extract.py $x'_length.faa' ./Filter/IPSCAN.tsv "PF00005" | sed 's/:/_/g' > $x'_length_NBD.faa'
grep -A 1 -f ./Filter/preliminary/Quality_threshold_species.txt $x'_length_NBD.faa' | sed '/--/d' > $x'_length_NBD_species.faa'
cat $x'_length_NBD_species.faa' | sed 's/\:/_/g' | sed 's/\-/_/g' | sed 's/\./_/g' > $x'_length_NBD_quality_clean.faa'

y=$x'_length_NBD_quality_clean.faa'
mafft --thread $THREADS $y > $y.aln
/home/pioannidis/Programs/trimAl/source/trimal -in $y'.aln' -out $y'.aln.trimm'
/data2/shane/Applications/custom/fasta_2_phylip.sh $y'.aln.trimm' > $y'.aln.trimm.phy'
Rscript /data2/shane/Applications/custom/Phylip_duplicate_NEW.R $y'.aln.trimm.phy' > $y'_clean.phy'
/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $y'_clean.phy' -n 'Unsorted.nwk'

Rscript ABC_include_unsorted.R