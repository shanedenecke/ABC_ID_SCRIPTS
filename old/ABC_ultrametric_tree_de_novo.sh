#!/usr/bin/env bash
H=/mnt/disk/shane/Transporter_ID/ABC_id
PHYLO=$H/ABC_REF/Input_files/Phylo_list.txt
SPEC=$H/ABC_REF/Input_files/target_species.tsv
QUAL_THRESH=.2
THREADS=14

cd $H

mkdir ./CAFE/Ultrametric_tree
mkdir ./CAFE/species_lists


## Create ultrametric Trees
for i in ./ABC_REF/CAFE/*.txt
do
  b=$(echo $(basename $i) | sed 's/.+CAFE.//g' | sed 's/_final_species.txt//g')
  mkdir ./CAFE/$b
  mkdir ./CAFE/$b/temp_proteome
  
  cat $i | while read i; do cp ./proteomes/$i* ./CAFE/$b/temp_proteome/; done
  
  ~/Applications/Orthofinder/orthofinder -f ./CAFE/temp_proteome -o ./CAFE/$b/orthofinder
  
  cat ./CAFE/temp_proteome/*.faa > ./CAFE/$b/total_proteome.faa
  python3 ~/Applications/Custom_Applications/Orthofind_parse.py -outdir ./CAFE/$b/one_to_one -inputdir ./CAFE/$b/orthofinder -total_fasta ./CAFE/$b/total_proteome.faa
  
  otodir=$(readlink -f ./CAFE/$b/one_to_one)
  
  ### perform alignments for all one to ones 
  for x in  $otodir/*.faa
  do
    cat $x | sed 's/\./_/g' | mafft --quiet --thread $THREADS - > $x'.aln'
    ~/Applications/trimal/source/trimal -automated1 -phylip_paml -in $x.aln -out $x.phy
  done
  
  
  #### merge all phylip files
  Rscript ~/Applications/Custom_Applications/Phylip_merge.R $otodir > $otodir'/Full_species.phy'
  
  #### make trees 
  if [ $b = 'Lepidopteran' ]; then
    echo 'lep'
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $otodir/Full_species.phy -n $b.tre -w ./CAFE/$b -o 286706_0
  elif [ $b = "Hemipteran" ]; then
    sed -i 's/J/A/g' $otodir/Full_species.phy
    sed -i 's/\./A/g' $otodir/Full_species.phy
    echo "Hemi"
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $otodir/Full_species.phy -n $b.tre -w ./CAFE/$b -o 7227_0
  elif [ $b = 'Diptera' ]; then
    ~/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $THREADS -m PROTGAMMAAUTO -s $otodir/Full_species.phy -n $b.tre -w ./CAFE/$b
  fi
  
  cp ./CAFE/$b   ./CAFE/clean_raxml_trees/
done



#raxml-ng --msa ABCF_.phy --prefix NG_test/test --all --model LG+G8+F --redo --threads 4 --bs-metric fbp,tbe --bs-trees autoMRE
#raxml-ng --all --msa prim.phy --model GTR+G --prefix T15  --threads 2 --bs-metric fbp,tbe
