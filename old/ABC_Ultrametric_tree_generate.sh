#!/usr/bin/env bash
H=~/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H


mkdir ./CAFE/Ultrametric_tree
mkdir ./CAFE/species_lists



for i in ./ABC_REF/CAFE/*.txt; do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_preliminary_species.txt//g')
  cp $i ./CAFE/species_lists/$b'_preliminary_species.txt'
  #grep -f ./Filter/preliminary/Quality_threshold_taxid.txt ./CAFE/species_lists/$b'_preliminary_species.txt' > ./CAFE/species_lists/$b'_final_species.txt'
  cat ./CAFE/species_lists/$b'_preliminary_species.txt' > ./CAFE/species_lists/$b'_final_species.txt'
done


## Create ultrametric Trees
for i in ./CAFE/species_lists/*final*
do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_final_species.txt//g')
  
  if [ $b = "Arthropod" ]; then
    ~/Applications/Custom_Applications/one_to_one_ID_exec.py -node "Arthropod" -taxid $i -output seq
  else
   ~/Applications/Custom_Applications/one_to_one_ID_exec.py -node "Arthropod" -taxid $i -output seq
  fi
  
  ### rename identified foler of sequences 
  fulltemp=$(readlink -f './CAFE/Ultrametric_tree/'$b'_og_sequences')
  mv og_sequences $fulltemp
  
  
  ### perform alignments for all one to ones 
  for x in  $fulltemp/*.faa
  do
    cat $x | sed 's/\./_/g' | mafft --quiet --thread $THREADS - > $x'.aln'
    ~/Applications/trimAl/source/trimal -automated1 -in $x'.aln' -out $x'.aln.trimm'
    ~/Applications/Custom_Applications/fasta_2_phylip.sh $x'.aln.trimm' | sed '1d' > $x'.aln.trimm.phy'
  done
  
  #### merge all phylip files
  Rscript ~/Applications/Custom_Applications/Phylip_merge.R $fulltemp > $fulltemp'/Full_species.phy'
  
  #### make trees 
  if [ $b = 'Lepidopteran' ]; then
    echo 'lep'
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp -o 286706_0
  elif [ $b = "Hemipteran" ]; then
    sed -i 's/J/A/g' $fulltemp/Full_species.phy
    sed -i 's/\./A/g' $fulltemp/Full_species.phy
    echo "Hemi"
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp -o 7227_0
  elif [ $b = 'Diptera' ]; then
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b.tre -w $fulltemp
  fi
  
  cp $fulltemp'/RAxML_bipartitions.'$b'.tre' ./CAFE/clean_raxml_trees/
done
