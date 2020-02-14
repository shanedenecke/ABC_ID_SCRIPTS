#!/usr/bin/env bash
cd $H



#### get NBD of identfied ABCs
cat ./preliminary_ABC/proteomes/* >> ./ABC_phylo/ABC_total.faa
/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i ./ABC_phylo/ABC_total.faa -o ./ABC_phylo/IPSCAN.tsv -f TSV
/data2/shane/Applications/custom/ip_domain_extract.py ./ABC_phylo/ABC_total.faa ./ABC_phylo/IPSCAN.tsv "PF00005" > ./ABC_phylo/ABC_total_NBD.faa
cat ./model_database/model_NBDs/Only_ABCs_NBD.faa >> ./ABC_phylo/ABC_total_NBD.faa


### filter for species 
cat ./ABC_REF/Input_files/Phylo_list.txt | while read i; do
  grep -A 1 $i ./ABC_phylo/ABC_total_NBD.faa >> ./ABC_phylo/ABC_total_NBD_species_filter.faa
done

### sort into families 
cat ./ABC_REF/Input_files/ABC_families.txt | while read i; do 
  grep -A 1 $i ./ABC_phylo/ABC_total_NBD_species_filter.faa | sed '/--/d' > ./ABC_phylo/ABC_phylo_$i.faa
  mafft --thread $THREADS ./ABC_phylo/ABC_phylo_$i.faa > ./ABC_phylo/ABC_phylo_$i.faa.aln
  /data2/shane/Transporter_ID/ABC_id/ABC_ID_SCRIPTS/general_scripts/trimAl/source/trimal -in ./ABC_phylo/ABC_phylo_$i.faa.aln -out ./ABC_phylo/ABC_phylo_$i.faa.aln.trimm
  #/data2/shane/Transporter_ID/ABC_id/ABC_ID_SCRIPTS/general_scripts/fasta_2_phylip.sh ./ABC_phylo/ABC_phylo_$i.faa.aln.trimm ./ABC_phylo/ABC_phylo_$i.faa.aln.trimm.phy
 /data2/shane/Applications/custom/fasta_2_phylip.sh ./ABC_phylo/ABC_phylo_$i.faa.aln.trimm > ./ABC_phylo/ABC_phylo_$i.faa.aln.trimm.phy
done

### move into subdirectory and remove duplicates
mkdir ./ABC_phylo/trees/
mv ./ABC_phylo/*.phy ./ABC_phylo/trees/
Rscript ./ABC_ID_SCRIPTS/Phylip_duplicate.R $H/ABC_phylo/trees/

### make trees 
for i in ./ABC_phylo/trees/*.phy
do
  b=$(echo $(basename $i) | sed 's/_.faa.aln.trimm.phy//g') 
  raxfile=$i
  raxdir=$H/ABC_phylo/trees/
  $H/ABC_ID_SCRIPTS/general_scripts/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T $THREADS -m PROTGAMMAAUTO -s $raxfile -n $b'.tre' -w $raxdir 
done 

mkdir ./ABC_phylo/clean_trees
mv ./ABC_phylo/*bipartitions.* ./ABC_phylo/clean_trees/
