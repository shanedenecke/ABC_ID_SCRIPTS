#!/usr/bin/env R

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))

#### set up directories and args
#setwd('~/Transporter_ID/ABC_id')
dir.create('./Filter/preliminary',showWarnings = F)
dir.create('./Filter/Final_transporters',showWarnings = F)
dir.create('./Filter/Final_transporters/dicts',showWarnings = F)
dir.create('./Filter/Final_transporters/proteomes',showWarnings = F)
dir.create('./Filter/Total_combined',showWarnings = F)


### Import arguments
args = commandArgs(trailingOnly=TRUE)
#args[1]=.3
thresh=as.numeric(args[1])

### Import key and metadata
key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted'),domains=c(2,2,1,2,1,2,2,1,1,1))
metadata=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code)

################ FUNCTIONS

shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything())
  return(b)
}

domain.annot=function(domain.table){
  domain.table$species=gsub('(^.+)__.+__.+$','\\1',domain.table$name)
  domain.table$code=gsub('^.+__.+__(.+$)','\\1',domain.table$name)
  return(domain.table)
}

count.fams=function(x){
	y=x %>% group_by(species,fam) %>% summarize(count=length(code)) %>% arrange(species,fam) %>% data.table()
	return(y)
}

### Parse IPscan
pfam.scan=fread('./Filter/HMM_PF00005_output_clean.tsv',select=c(1,3,20,21)) %>% rename(name=V1,len=V3,start=V20,end=V21)
pfam.scan$fam=gsub('^.+__(.+)__.+$','\\1',pfam.scan$name)
pfam.table=pfam.scan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% filter(fam!="ABC_Unsorted_") %>% data.table()


#ipscan=fread('/data2/shane/Transporter_ID/ABC_id/Filter/old_IP/IP_combined.tsv',sep='\t',fill=T)[V5=='PF00005'] %>% select(V1,V3,V7,V8) %>% rename(name=V1,len=V3,start=V7,end=V8)
#ipscan=fread('./Filter/IPSCAN.tsv',sep='\t',fill=T)[V5=='PF00005'] %>% select(V1,V3,V7,V8) %>% rename(name=V1,len=V3,start=V7,end=V8)
#ipscan$fam=gsub('^.+__(.+)__.+$','\\1',ipscan$name)
#iptable=ipscan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% filter(fam!="ABC_Unsorted_") %>% data.table()



total.fasta=seqinr::read.fasta('./Filter/ABC_preliminary_total.faa',set.attributes = F,as.string = T,forceDNAtolower = F)

short.list=list()
good.list=list()
long.list=list()
for(i in unique(pfam.table$fam)){
  sub=pfam.table[fam==i]
  mins=key[family==i]$domains
  short.list[[i]]=sub[domains<mins]
  good.list[[i]]=sub[domains==mins]
  long.list[[i]]=sub[domains>mins]
}
too.short=rbindlist(short.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()
too.long=rbindlist(long.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()
just.right=rbindlist(good.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()

### get counts
short.count.sp=too.short %>% count.fams() %>% group_by(species) %>% summarize(short.count=sum(count)) %>% data.table()
long.count.sp=too.long %>% count.fams() %>% group_by(species) %>% summarize(long.count=sum(count)) %>% data.table()
good.count.sp=just.right %>% count.fams() %>% group_by(species) %>% summarize(good.count=sum(count)) %>% data.table()


##### Quality cutoff
m=merge(short.count.sp,long.count.sp,by='species',all=T) %>% 
  merge(good.count.sp,by='species',all=T) %>%  replace(is.na(.), 0) %>% 
  mutate(frac_bad=(short.count+long.count)/(short.count+long.count+good.count)) %>%
  data.table() 
good.species=m[frac_bad<thresh]$species
quality.table=m %>% rename(abbreviation=species) %>% merge(metadata,by='abbreviation') 
good.taxid=quality.table[frac_bad<thresh]$taxid_code


fwrite(merge(quality.table,metadata,by='abbreviation'),'./Filter/preliminary/Quality_table.tsv',sep='\t')
writeLines(good.species,'./Filter/preliminary/Quality_threshold_species.txt')
writeLines(good.taxid,'./Filter/preliminary/Quality_threshold_taxid.txt')




########### Final outputs

fwrite(too.short,'./Filter/preliminary/pelim_short.csv')
fwrite(too.long,'./Filter/preliminary/pelim_long.csv')
fwrite(just.right,'./Filter/preliminary/pelim_good.csv')



