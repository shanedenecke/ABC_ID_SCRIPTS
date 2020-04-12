#!/usr/bin/env bash

mkdir phylo
mkdir ./phylo/clean_trees

#### get NBD of identfied ABCs
grep -f ./ABC_REF/Input_files/Phylo_list.txt -A 1 ./Final_outputs/Combined_files/Total_ABC.faa | sed '/--/d' > ./phylo/ABC_total.faa
python3 ~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABC_total.faa  -table ./Filter/HMM_PF00005_output.tsv > ./phylo/ABC_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABC_total.faa
sed -i 's/\./_/g' ./phylo/ABC_total.faa



##### ONly ABCH

mkdir ./phylo/ABCH
########## Make ABCH tree specifically
grep -E -A 1 'DroMel|NezVir|BemTab|MyzPer|NilLug|TetUrt' ./Final_outputs/Combined_files/Total_ABC.faa | grep -A 1 'ABCH' | sed '/--/d' > ./phylo/ABCH/ABCH_total.faa
#python3 ~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABCH/ABCH_total.faa  -table ./Filter/HMM_PF00005_output.tsv > ./phylo/ABCH/ABCH_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABCH/ABCH_total.faa
sed -i 's/\./_/g' ./phylo/ABCH/ABCH_total.faa

x=./phylo/ABCH/ABCH_total
cat ./ABC_REF/phylo_premade/EscCol_ARTP_outgroup_full.faa >> $x.faa
mafft --quiet --thread $THREADS $x.faa > $x.aln
~/Applications/trimal/source/trimal -strict -phylip_paml -in $x.aln -out $x.phy


~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $x.phy -n ABCH.nwk -w $H/phylo/ABCH/ -o Ecoli_ARTP_OUT



### sort into families
x=./phylo/ABC_total
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 '__'$i ./phylo/ABC_total.faa | sed '/--/d' > ./phylo/$i.faa
  cat ./ABC_REF/phylo_premade/EscCol_ARTP_outgroup_full.faa >> ./phylo/$i.faa
  mafft --quiet --thread $THREADS ./phylo/$i.faa > ./phylo/$i.aln
  ~/Applications/trimal/source/trimal -strict -phylip_paml -in ./phylo/$i.aln -out ./phylo/$i.phy 
done


### make trees All 
for i in ./phylo/*.phy
do
  b=$(echo $(basename $i) | sed 's/_.phy//g')
  out=$(grep 'OUT' './phylo/'$b'_.faa' | sed 's/>//g')
  #raxfile=$H'/phylo/'$(basename $i)
  raxdir=$H/phylo/
   ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $i -n $b'.nwk' -w $raxdir -o $out
done 


