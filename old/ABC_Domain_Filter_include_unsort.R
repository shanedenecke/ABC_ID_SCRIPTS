#!/usr/bin/env R

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))
shhh(library(ggtree))
shhh(library(treeio))
shhh(library(ggplot2))


### Import key and metadata
#setwd('/home/sdenecke/Transporter_ID/ABC_id')
dir.create('./Filter/Final_transporters',showWarnings = F)
dir.create('./Filter/Final_transporters/dicts',showWarnings = F)
dir.create('./Filter/Final_transporters/proteomes',showWarnings = F)
dir.create('./Filter/Total_combined',showWarnings = F)


key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted'),domains=c(2,2,1,2,1,2,2,1,1,1))
metadata=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code)
good.species=readLines('./Filter/preliminary/Quality_threshold_species.txt')
total.fasta=seqinr::read.fasta('./Filter/ABC_preliminary_total.faa',set.attributes = F,as.string = T,forceDNAtolower = F)
################ FUNCTIONS

domain.annot=function(domain.table){
  domain.table$species=gsub('(^.+)__.+__.+$','\\1',domain.table$name)
  domain.table$code=gsub('^.+__.+__(.+$)','\\1',domain.table$name)
  return(domain.table)
}

count.fams=function(x){
  y=x %>% group_by(species,fam) %>% summarize(count=length(code)) %>% arrange(species,fam) %>% data.table()
  return(y)
}

### Parse PfamScan
pfam.scan=fread('./Filter/HMM_PF00005_output_clean.tsv',select=c(1,3,20,21)) %>% rename(name=V1,len=V3,start=V20,end=V21)
pfam.scan$fam=gsub('^.+__(.+)__.+$','\\1',pfam.scan$name)
pfam.table=pfam.scan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% filter(fam!="ABC_Unsorted_") %>% data.table()



good.list=list()
for(i in unique(pfam.table$fam)){
  sub=pfam.table[fam==i]
  mins=key[family==i]$domains
  good.list[[i]]=sub[domains==mins]
}
good=rbindlist(good.list) %>% domain.annot() %>% filter(species %in% good.species) %>% data.table()



l=good %>% arrange(species) %>% group_split(species)
sps=sort(good$species) %>% unique()
for(i in 1:length(l)){
  fwrite(l[[i]],file=paste0('./Filter/Final_transporters/dicts/',sps[i],'_final_ABC_dict.csv'))
  nam=l[[i]]$name
  sub.fa=total.fasta[names(total.fasta) %in% nam]
  write.fasta(sub.fa,names=names(sub.fa),file.out = paste0('./Filter/Final_transporters/proteomes/',sps[i],'_final.faa'))
}


######################### TOTAL OUTPUTS


#### produce counts files
good.sum=good %>% count.fams() %>% spread(key=species,value=count) %>%
  data.table()
good.sum[is.na(good.sum)]=0
good.sum=data.table(fam=good.sum$fam,apply(select(good.sum,2:length(colnames(good.sum))),2,as.numeric)) 
fwrite(good.sum,'./Filter/Total_combined/Full_transporter_counts.csv') 


#### write fasta
abc.fasta=total.fasta[names(total.fasta) %in% good$name]
write.fasta(abc.fasta,names=names(abc.fasta),file.out = './Filter/Total_combined/Total_ABC.faa')

### write dictionary
fwrite(good,'./Filter/Total_combined/Total_ABC_dict.csv')
