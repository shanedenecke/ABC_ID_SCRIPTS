cd $H

mkdir memetemp
cat ./GENERAL_REFERENCE/Model_ABC_sets/ABC_proteins/* > ./memetemp/All_models.faa
meme All_models.faa -mod zoops -nmotifs 3 -o ABC_meme_motif -protein -minw 4 -maxw 8


touch ./ABC_meme_motif/meme_clean.txt
echo "MEME version 4ALPHABET= ACDEFGHIKLMNPQRSTVWY" >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 3 "Background letter frequencies" >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep "MEME-1 position-specific probability matrix" | perl -pe 's/Motif ([A-Z]+) .+$/MOTIF $1/g' >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 2 "MEME-1 position-specific probability matrix" | tail -n 1 >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 10 "MEME-1 position-specific probability matrix" | tail -n 8 >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep "MEME-2 position-specific probability matrix" | perl -pe 's/Motif ([A-Z]+) .+$/MOTIF $1/g' >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 2 "MEME-2 position-specific probability matrix" | tail -n 1 >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 10 "MEME-2 position-specific probability matrix" | tail -n 8 >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep "MEME-3 position-specific probability matrix" | perl -pe 's/Motif ([A-Z]+) .+$/MOTIF $1/g' >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 2 "MEME-3 position-specific probability matrix" | tail -n 1 >> ./ABC_meme_motif/meme_clean.txt
cat ./ABC_meme_motif/meme.txt | grep -A 10 "MEME-3 position-specific probability matrix" | tail -n 8 >> ./ABC_meme_motif/meme_clean.txt
sed -i 's/\tMOTIF/MOTIF/g' ./ABC_meme_motif/meme_clean.txt 


mast -hit_list ./ABC_meme_motif/meme_clean.txt ../ABC_search/HelArm/Preliminary_ABCs.faa -notext -nohtml -nostatus | tail -n +3 | head -n -1 | sed -E 's/\s+/___/g' | sed 's/___/\t/g' | cut -f 1,3
