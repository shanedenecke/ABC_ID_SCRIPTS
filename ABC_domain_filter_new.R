#!/usr/bin/env R
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))

#setwd('~/Transporter_ID/ABC_id/')
dir.create('./Filter',showWarnings = F)
dir.create('./Final_outputs',showWarnings = F)
dir.create('./Final_outputs/Combined_files',showWarnings = F)
dir.create('./Final_outputs/proteomes',showWarnings = F)
dir.create('./Final_outputs/dicts',showWarnings = F)


system('cat ./preliminary_ABC/proteomes/* > ./Filter/ABC_preliminary_total.faa')
system('cat ./Db_build_temp/Only_ABCs.faa* >> ./Filter/ABC_preliminary_total.faa')
system('hmmsearch --domtblout ./Filter/HMM_PF00005_output.tsv ./ABC_REF/Model_ABC_sets/ABC_tran.hmm ./Filter/ABC_preliminary_total.faa > ./Filter/hmm_junk.txt')
system('cat ./Filter/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > ./Filter/HMM_PF00005_output_clean.tsv')


### Import arguments
args = commandArgs(trailingOnly=TRUE)
#args[1]=.3
thresh=as.numeric(args[1])


### Import key and metadata
key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted'),
               mins=c(1,2,1,2,1,2,2,1,1,1),maxes=c(c(2,2,1,2,1,2,2,1,1,2)))
metadata=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code)
total.fasta=seqinr::read.fasta('./Filter/ABC_preliminary_total.faa',set.attributes = F,as.string = T,forceDNAtolower = F)
n50=fread('./ABC_REF/species_metadata/i5k_n50.csv')

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
  dom.max=key[family==i]$maxes
  dom.min=key[family==i]$min
  short.list[[i]]=sub[dom.min<mins & len<150]
  good.list[[i]]=sub[domains==mins & len>150]
  #long.list[[i]]=sub[domains>mins]
}
too.short=rbindlist(short.list) %>% domain.annot() %>% data.table()
too.long=rbindlist(long.list) %>% domain.annot() %>% data.table()
just.right=rbindlist(good.list) %>% domain.annot() %>% data.table()


### get counts
short.count.sp=too.short %>% count.fams() %>% group_by(species) %>% summarize(short.count=sum(count)) %>% data.table()
#long.count.sp=too.long %>% count.fams() %>% group_by(species) %>% summarize(long.count=sum(count)) %>% data.table()
good.count.sp=just.right %>% count.fams() %>% group_by(species) %>% summarize(good.count=sum(count)) %>% data.table()


##### Quality cutoff
quality.table=merge(short.count.sp,good.count,by='species',all=T) %>% replace(is.na(.), 0) %>% 
  mutate(frac_bad=(short.count/(short.count+good.count))) %>%
  rename(abbreviation=species) %>% merge(metadata,by='abbreviation') %>%
  merge(n50,by='Species_name',all.x=T) %>%
  data.table() 

#q2=quality.table[frac_bad<.4 & (N50>50000 | is.na(N50))]
q2=quality.table[frac_bad<thresh]

good.species=q2[frac_bad<thresh]$abbreviation
good.taxid=q2[frac_bad<thresh]$taxid_code


fwrite(quality.table,'./Filter/Quality_table.tsv',sep='\t')
writeLines(good.species,'./Filter/Quality_threshold_species.txt')
writeLines(good.taxid,'./Filter/Quality_threshold_taxid.txt')


########### Output good
good.l=just.right %>% arrange(species) %>% group_split(species)
for(i in 1:length(good.l)){
  sp=good.l[[i]]$species[1]
  fwrite(good.l[[i]],file=paste0('./Final_outputs/dicts/',sp,'_final_ABC_dict.csv'))
  nam=good.l[[i]]$name
  sub.fa=total.fasta[names(total.fasta) %in% nam]
  write.fasta(sub.fa,names=names(sub.fa),file.out = paste0('./Final_outputs/proteomes/',sp,'_final.faa'))
}


######################### TOTAL OUTPUTS


#### produce counts files
good.sum=just.right %>% count.fams() %>% spread(key=species,value=count) %>%
  data.table()
good.sum[is.na(good.sum)]=0
good.sum=data.table(fam=good.sum$fam,apply(select(good.sum,2:length(colnames(good.sum))),2,as.numeric)) 
fwrite(good.sum,'./Final_outputs/Combined_files/ABC_transporter_counts.csv') 


#### write fasta
abc.fasta=total.fasta[names(total.fasta) %in% just.right$name]
write.fasta(abc.fasta,names=names(abc.fasta),file.out = './Final_outputs/Combined_files/Total_ABC.faa')

### write dictionary
fwrite(just.right,'./Final_outputs/Combined_files/Total_ABC_dict.csv')

