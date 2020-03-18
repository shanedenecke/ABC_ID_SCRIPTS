cd $H




## create CAFE script files
mkdir ./CAFE/scripts
mkdir ./CAFE/logfiles
mkdir ./CAFE/outputs

for i in ./CAFE/CAFE_tables/*.tsv
do
  b=$(echo $(basename $i) | sed 's/_ABC_CAFE_table.tsv//g')
  
  ##### Make CAFE Script
  rm ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  touch ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  echo '#!cafe' >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  echo 'load -i '$H'/CAFE/CAFE_tables/'$b'_ABC_CAFE_table.tsv -t '$THREADS' -l '$H'/CAFE/logfiles/'$b'_ABC_logfile.txt -p 0.05' >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  a=$(cat ./CAFE/clean_raxml_trees/$b'_tree_ultrametric.tre')
  echo 'tree '$a >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  l=$(cat ./CAFE/clean_raxml_trees/$b'_tree_lambda.txt')
  #echo 'lambda -l '$lambda' -t '$l >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  echo 'lambda -s -t '$l >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  echo 'report '$H'/CAFE/outputs/'$b'_ABC_cafe_output' >> ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  
  ## run cafe
  ~/Applications/CAFE/release/cafe ./CAFE/scripts/$b'_ABC_Cafe_Script.cafe'
  python2  ./ABC_ID_SCRIPTS/Fulton_python_scripts/cafetutorial_report_analysis.py -i ./CAFE/outputs/$b'_ABC_cafe_output.cafe' -o ./CAFE/outputs/$b'_ABC_summary.txt'


  t=$(grep 'Tree:' ./CAFE/outputs/$b'_ABC_cafe_output.cafe' | sed 's/Tree://g')
  d=$(grep "IDs of nodes:" ./CAFE/outputs/$b'_ABC_cafe_output.cafe' | sed 's/# IDs of nodes://g')
  
  #python2 ./ABC_ID_SCRIPTS/Fulton_python_scripts/cafetutorial_draw_tree.py -i ./CAFE/outputs/$b'_ABC_summary.txt_node.txt' \
  #-t $t \
  #-d $d \
  #-y Expansions -o ./CAFE/outputs/$b'_Expansions_tree.png'
done