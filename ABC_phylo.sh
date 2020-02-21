#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/ABC_id'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=10

cd $H


mkdir phylo

#### get NBD of identfied ABCs
rm ./phylo/ABC_total.faa
cat ./Filter/Full_transporters/proteomes/* >> ./phylo/ABC_total.faa
/data2/shane/Applications/custom/ip_domain_extract.py ./phylo/ABC_total.faa ./Filter/IPSCAN.tsv "PF00005" > ./phylo/ABC_total_NBD.faa
#cat ./model_database/model_NBDs/Only_ABCs_NBD.faa >> ./phylo/ABC_total_NBD.faa


### filter for species 
cat ./ABC_REF/Input_files/Phylo_list.txt | while read i; do
  grep -A 1 $i ./phylo/ABC_total_NBD.faa >> ./phylo/ABC_total_NBD_species_filter.faa
done
sed -i 's/-/_/g' ./phylo/ABC_total_NBD_species_filter.faa

### sort into families 
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 '__'$i ./phylo/ABC_total_NBD_species_filter.faa | sed '/--/d' > ./phylo/phylo_$i.faa
  grep -A 1 $i ./ABC_REF/Model_ABC_sets/ABC_outgroups.faa | sed '/--/d' >> ./phylo/phylo_$i.faa
  mafft --thread $THREADS ./phylo/phylo_$i.faa > ./phylo/phylo_$i.faa.aln
  /data2/shane/Transporter_ID/ABC_id/ABC_ID_SCRIPTS/general_scripts/trimAl/source/trimal -in ./phylo/phylo_$i.faa.aln -out ./phylo/phylo_$i.faa.aln.trimm
  #/data2/shane/Transporter_ID/ABC_id/ABC_ID_SCRIPTS/general_scripts/fasta_2_phylip.sh ./phylo/phylo_$i.faa.aln.trimm ./phylo/phylo_$i.faa.aln.trimm.phy
 /data2/shane/Applications/custom/fasta_2_phylip.sh ./phylo/phylo_$i.faa.aln.trimm > ./phylo/phylo_$i.faa.aln.trimm.phy
done

### move into subdirectory and remove duplicates
mkdir ./phylo/trees/
mv ./phylo/*.phy ./phylo/trees/
Rscript ./ABC_ID_SCRIPTS/Phylip_duplicate.R $H/phylo/trees/

mkdir ./phylo/clean_trees
### make trees 
for i in ./phylo/trees/*.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy//g')
  mkdir ./phylo/clean_trees/$b
  out=$(grep 'outgroup' './phylo/'$b'_.faa' | sed 's/>//g')
  raxfile=$i
  raxdir=$H/phylo/clean_trees/$b
  $H/ABC_ID_SCRIPTS/general_scripts/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir -o $out
done 

