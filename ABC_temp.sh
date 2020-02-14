#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/ABC_id'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SLC_FAM=$H/ABC_REF/Input_files/ABC_Families.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
THREADS=12

cd $H


### phylogeny
### make trees 
for i in ./ABC_phylo/trees/*.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy//g') 
  raxfile=$i
  raxdir=$H/ABC_phylo/trees/
  $H/ABC_ID_SCRIPTS/general_scripts/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 200 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
done 

mkdir ./ABC_phylo/clean_trees
mv ./ABC_phylo/*bipartitions.* ./ABC_phylo/clean_trees/