#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/ABC_id'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H

#### get NBD of identfied ABCs
rm ./phylo/ABC_total.faa
grep -f ./ABC_REF/Input_files/Phylo_list.txt -A 1 ./Filter/Total_combined/Total_ABC.faa > ./phylo/ABC_total.faa
/data2/shane/Applications/custom/ip_domain_extract.py ./phylo/ABC_total.faa ./Filter/IPSCAN.tsv "PF00005" > ./phylo/ABC_total_NBD.faa

sed -i 's/-/_/g' ./phylo/ABC_total_NBD.faa
sed -i 's/\./_/g' ./phylo/ABC_total_NBD.faa


### sort into families
x=./phylo/ABC_total_NBD
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 '__'$i $x.faa | sed '/--/d' > $x$i.faa
  grep -A 1 $i ./ABC_REF/Model_ABC_sets/ABC_outgroups.faa | sed '/--/d' >> $x$i.faa
  mafft --quiet --thread $THREADS $x$i.faa > $x$i.faa.aln
  /data2/shane/Transporter_ID/ABC_id/ABC_ID_SCRIPTS/general_scripts/trimAl/source/trimal -in $x$i.faa.aln -out $x$i.faa.aln.trimm
 /data2/shane/Applications/custom/fasta_2_phylip.sh  $x$i.faa.aln.trimm > $x$i.faa.aln.trimm.phy
 Rscript /data2/shane/Applications/custom/Phylip_duplicate.R $x$i.faa.aln.trimm.phy > $x$i.faa.aln.trimm.phy.phy
done


### make trees 
for i in ./phylo/*.phy.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy.phy//g')
  #out=$(grep 'outgroup' './phylo/'$b'_.faa' | sed 's/>//g')
  raxfile=$H'/phylo/'$(basename $i)
  raxdir=$H/phylo/clean_trees/
  ./ABC_ID_SCRIPTS/general_scripts/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.nwk' -w $raxdir
done 



