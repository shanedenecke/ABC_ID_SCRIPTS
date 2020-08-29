#!/usr/bin/env bash

H=/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
SPEC=$H/GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv
BUSCO_THRESH=75
THREADS=6

### 2) Set Home working Directory
cd $H 



mkdir phylo
mkdir ./phylo/clean_trees

#### get NBD of identfied ABCs
rm ./phylo/Phylo_list.txt
for i in SpoFru BomMor PapPol DanPle HalHal; do echo $i >> ./phylo/Phylo_list.txt; done

grep -f ./phylo/Phylo_list.txt -A 1 ./Final_outputs/combined_files/All_ABCs.faa | sed '/--/d' > ./phylo/ABC_total.faa
python3 ~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABC_total.faa  -table ./model_database/Db_build_temp/PF0005_HMM_output.tsv > ./phylo/ABC_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABC_total.faa
sed -i 's/\./_/g' ./phylo/ABC_total.faa


### sort into families
x=./phylo/ABC_total
for i in ABCBF; do 
  grep -A 1 '__'$i ./phylo/ABC_total.faa | sed '/--/d' > ./phylo/$i.faa
  #cat ./ABC_REF/phylo_premade/EscCol_ARTP_outgroup_full.faa >> ./phylo/$i.faa
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
  ~/Applications/raxml-ng --all --msa $i --prefix ./phylo/$b'_tree' --threads 6  --bs-trees autoMRE{1000} --model LG+G8+F --redo 
done 



##### ONly ABCH

mkdir ./phylo/ABCH
########## Make ABCH tree specifically
grep -E -A 1 'NezVir|DakVit|MyzPer|NilLug|DapPul|TetUrt' ./Final_outputs/combined_files/All_ABCs.faa | grep -A 1 'ABCH' | sed '/--/d' > ./phylo/ABCH/ABCH_total.faa
#python3 ~/Applications/Custom_Applications/hmmsearch_pfam_domain_parse.py -fasta ./phylo/ABCH/ABCH_total.faa  -table ./Filter/HMM_PF00005_output.tsv > ./phylo/ABCH/ABCH_total_NBD.faa
sed -i 's/-/_/g' ./phylo/ABCH/ABCH_total.faa
sed -i 's/\./_/g' ./phylo/ABCH/ABCH_total.faa

x=./phylo/ABCH/ABCH_total
#cat ./ABC_REF/phylo_premade/EscCol_ARTP_outgroup_full.faa >> $x.faa
#mafft --quiet --thread $THREADS $x.faa > $x.aln
#~/Applications/trimal/source/trimal -strict -phylip_paml -in $x.aln -out $x.phy

#~/Applications/raxml-ng --all --msa $x.phy --prefix $x'_tree' --threads 6  --bs-trees autoMRE{1000} --model LG+G8+F --redo 





