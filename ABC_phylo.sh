#!/usr/bin/env bash

#### get NBD of identfied ABCs
grep -f ./ABC_REF/Input_files/Phylo_list.txt -A 1 ./Final_outputs/Combined_files/All_species_ABCS.faa | sed '/--/d' > ./phylo/ABC_total.faa
python3 ~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABC_total.faa  -table ./Filter/HMM_PF00005_output.tsv > ./phylo/ABC_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABC_total_NBD.faa
sed -i 's/\./_/g' ./phylo/ABC_total_NBD.faa


### sort into families
x=./phylo/ABC_total_NBD
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 '__'$i $x.faa | sed '/--/d' > $x$i.faa
  #grep -A 1 $i ./ABC_REF/Model_ABC_sets/ABC_outgroups.faa | sed '/--/d' >> $x$i.faa
  cat ./ABC_REF/phylo_premade/EscCol_ARTP_outgroup.faa >> $x$i.faa
  mafft --quiet --thread $THREADS $x$i.faa > $x$i.faa.aln
  ~/Applications/trimAl/source/trimal -automated1 -in $x$i.faa.aln -out $x$i.faa.aln.trimm 
 ~/Applications/Custom_Applications/fasta_2_phylip.sh  $x$i.faa.aln.trimm > $x$i.faa.aln.trimm.phy
 Rscript ~/Applications/Custom_Applications/Phylip_duplicate.R $x$i.faa.aln.trimm.phy > $x$i.faa.aln.trimm.phy.phy
done


### make trees 
for i in ./phylo/*.phy.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy.phy//g')
  out=$(grep 'OUT' './phylo/'$b'_.faa' | sed 's/>//g')
  #raxfile=$H'/phylo/'$(basename $i)
  raxdir=$H/phylo/clean_trees/
   ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $i -n $b'.nwk' -w $raxdir -o $out
done 


