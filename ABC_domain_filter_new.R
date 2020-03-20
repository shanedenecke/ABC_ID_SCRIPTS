#!/usr/bin/env R
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))


system('cat ./preliminary_ABC/proteomes/* > ./Filter/ABC_preliminary_total.faa')
system('cat ./Db_build_temp/Only_ABCs.faa* >> ./Filter/ABC_preliminary_total.faa')
system('hmmsearch --domtblout ./Filter/HMM_PF00005_output.tsv ./ABC_REF/Model_ABC_sets/ABC_tran.hmm ./Filter/ABC_preliminary_total.faa > ./Filter/hmm_junk.txt')
system('cat ./Filter/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > ./Filter/HMM_PF00005_output_clean.tsv')


### Import arguments
args = commandArgs(trailingOnly=TRUE)
#args[1]=.3
thresh=as.numeric(args[1])


### Import key and metadata
key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted'),domains=c(2,2,1,2,1,2,2,1,1,1))
metadata=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code)
total.fasta=seqinr::read.fasta('./Filter/ABC_preliminary_total.faa',set.attributes = F,as.string = T,forceDNAtolower = F)


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
pfam.table=pfam.scan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% data.table()

#### sort into categories
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
quality.table=merge(short.count.sp,long.count.sp,by='species',all=T) %>% 
  merge(good.count.sp,by='species',all=T) %>%  replace(is.na(.), 0) %>% 
  mutate(frac_bad=(short.count+long.count)/(short.count+long.count+good.count)) %>%
  rename(abbreviation=species) %>% merge(metadata,by='abbreviation')
  data.table() 

good.species=m[frac_bad<thresh]$species
good.taxid=quality.table[frac_bad<thresh]$taxid_code


fwrite(quality.table,'./Filter/Quality_table.tsv',sep='\t')
writeLines(good.species,'./Filter/Quality_threshold_species.txt')
writeLines(good.taxid,'./Filter/Quality_threshold_taxid.txt')


