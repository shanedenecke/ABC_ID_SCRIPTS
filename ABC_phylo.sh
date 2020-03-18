#!/usr/bin/env bash
H=~/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H

#### get NBD of identfied ABCs
grep -f ./ABC_REF/Input_files/Phylo_list.txt -A 1 ./Filter/Total_combined/Total_ABC.faa | sed '/--/d' > ./phylo/ABC_total.faa
#/data2/shane/Applications/custom/ip_domain_extract.py ./phylo/ABC_total.faa ./Filter/IPSCAN.tsv "PF00005" > ./phylo/ABC_total_NBD.faa
~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABC_total.faa  -table ./Filter/HMM_PF00005_output.tsv > ./phylo/ABC_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABC_total_NBD.faa
sed -i 's/\./_/g' ./phylo/ABC_total_NBD.faa


### sort into families
x=./phylo/ABC_total_NBD
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 '__'$i $x.faa | sed '/--/d' > $x$i.faa
  grep -A 1 $i ./ABC_REF/Model_ABC_sets/ABC_outgroups.faa | sed '/--/d' >> $x$i.faa
  mafft --quiet --thread $THREADS $x$i.faa > $x$i.faa.aln
  ~/Applications/trimAl/source/trimal -automated1 -in $x$i.faa.aln -out $x$i.faa.aln.trimm 
 ~/Applications/Custom_Applications/fasta_2_phylip.sh  $x$i.faa.aln.trimm > $x$i.faa.aln.trimm.phy
 Rscript ~/Applications/Custom_Applications/Phylip_duplicate.R $x$i.faa.aln.trimm.phy > $x$i.faa.aln.trimm.phy.phy
done


### make trees 
for i in ./phylo/*.phy.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy.phy//g')
  #out=$(grep 'outgroup' './phylo/'$b'_.faa' | sed 's/>//g')
  #raxfile=$H'/phylo/'$(basename $i)
  raxdir=$H/phylo/clean_trees/
   ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $i -n $b'.nwk' -w $raxdir
done 



