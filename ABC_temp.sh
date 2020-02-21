#!/usr/bin/env bash
H='/data2/shane/Transporter_ID/ABC_id'
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
THREADS=12

cd $H
mkdir CAFE
mkdir ./CAFE/orthoDB
mkdir ./CAFE/clean_raxml_trees

## Create ultrametric Trees
for i in ./ABC_REF/CAFE/*
do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_species_list.txt//g')
  
  if [ $b = "Arthropod" ] || [ $b = "Arachnid" ]; then
    /data2/shane/Applications/custom/OrthoDB/one_to_one_ID_exec.py -node "Metazoa" -taxid $i -output seq
  else
    /data2/shane/Applications/custom/OrthoDB/one_to_one_ID_exec.py -node "Arthropod" -taxid $i -output seq
  fi
  
  ### rename identified foler of sequences 
  temp='./CAFE/orthoDB/'$b'_og_sequences'
  mv og_sequences $temp
  
  
  ### perform alignments for all one to ones 
  for x in  $temp/*
  do
    mafft --thread $THREADS $x > $x'.aln'
    /data2/shane/Applications/trimAl/source/trimal -in $x'.aln' -out $x'.aln.trimm'
    /data2/shane/Applications/custom/fasta_2_phylip.sh $x'.aln.trimm' | sed '1d' > $x'.aln.trimm.phy'
  done
  #### merge all phylip files
  Rscript ./ABC_ID_SCRIPTS/Phylip_merge.R $temp
  
  fulltemp=$(readlink -f $temp)
  #### make trees 
  if [ $b = "Arthropod" ] || [ $b = "Arachnid" ]; then
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b'.tre' -w $raxdir -o 6239_0
  elif [ $b = 'Lepidopteran' ]; then
	#echo 'lep'       
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b'.tre' -w $raxdir -o 7029_0
  elif [ $b = "Hemipteran" ]; then
	sed -i 's/J/A/g' $fulltemp/Full_species.phy
	sed -i 's/\./A/g' $fulltemp/Full_species.phy
	/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $fulltemp/Full_species.phy -n $b'.tre' -w $raxdir -o 7227_0
 fi
 
 cp $raxdir/RAxML_bipartitions.$b.tre ./CAFE/clean_raxml_trees/$b'_raxml_clean.tre'
 
done
