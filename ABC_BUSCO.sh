cd $H

mkdir BUSCO
mkdir BUSCO/clean_summary

cd BUSCO
for i in ../proteomes/* ../ABC_REF/Model_ABC_sets/Raw_proteomes/*; do
  b=$(echo $(basename $i) | sed 's/_unigene.faa//g')
  mkdir ./BUSCO/$b
  
  if [ $b=='CaeEle' ] || [ $b=='HomSap' ]
  then
    python3 ~/Applications/busco/bin/busco -c $THREADS --quiet --config ~/Applications/busco/config/myconfig.ini -m proteins -i $i -o $b -f -l metazoa_odb10
  else
    python3 ~/Applications/busco/bin/busco -c $THREADS --quiet --config ~/Applications/busco/config/myconfig.ini -m proteins -i $i -o $b -f -l arthropoda_odb10
  fi 
  cp $b/*.txt ./clean_summary/$b'_clean_summary.txt'
  echo 'Finihsed with Species '$b
done

echo 'Finihsed with all species. Begginning to parse outputs'
python3 ~/Applications/Custom_Applications/BUSCO_parse.py -dir ./BUSCO/clean_summary/ > $H/BUSCO/BUSCO_final_summary.tsv


cd $H