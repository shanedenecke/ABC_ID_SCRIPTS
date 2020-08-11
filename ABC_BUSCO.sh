cd $H

mkdir BUSCO
mkdir BUSCO/clean_summary
mkdir BUSCO/junk

cd BUSCO
for i in ../proteomes/* ../GENERAL_REFERENCE/Model_ABC_sets/Raw_proteomes/*; do
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

## parse BUSCO outputs
~/Applications/Custom_Applications/BUSCO_parse.py -dir ./clean_summary/ -thresh $BUSCO_THRESH > $H/BUSCO/BUSCO_final_summary.tsv
~/Applications/Custom_Applications/BUSCO_parse.py -dir $H/BUSCO/clean_summary/ > $H/BUSCO/BUSCO_final_summary_unfiltered.tsv

ls $H/BUSCO/ | grep -E '^[[:alpha:]]{6}$' | while read i; do mv $H/BUSCO/$i* $H/BUSCO/junk/;done ### move all BUSCO junk to junk folder
#cat $H/BUSCO/BUSCO_final_summary.tsv | cut -f 6 | tail -n +2 | while read i; do cp $PROTEOMES/"$i"* $H/filtered_proteomes/; done ### copy all good proteomes to final folder
  
cd $H
