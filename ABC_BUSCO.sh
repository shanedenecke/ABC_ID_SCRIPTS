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
if grep -q 'TetUrt' $H/BUSCO/BUSCO_final_summary.tsv; then echo "Tetranychus already present"; else grep 'TetUrt' $H/BUSCO/BUSCO_final_summary_unfiltered.tsv >> $H/BUSCO/BUSCO_final_summary.tsv; fi

ls $H/BUSCO/ | grep -E '^[[:alpha:]]{6}$' | while read i; do mv $H/BUSCO/$i* $H/BUSCO/junk/;done ### move all BUSCO junk to junk folder

mkdir $H/proteomes/temp
mkdir $H/proteomes/low_busco_proteomes
cat $H/BUSCO/BUSCO_final_summary.tsv | cut -f 6 | tail -n +2 | while read i; do mv $H/proteomes/"$i"* $H/proteomes/temp; done ### copy all good proteomes to final folder

mv $H/proteomes/*.faa $H/proteomes/low_busco_proteoms/
mv $H/proteomes/temp/* $H/proteomes/
rm -r $H/proteomes/temp

cd $H
