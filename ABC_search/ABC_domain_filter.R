#!/usr/bin/env Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))


### Functions

domain.annot=function(domain.table){
  domain.table$abbreviation=gsub('(^.+)__.+__.+$','\\1',domain.table$name)
  domain.table$code=gsub('^.+__.+__(.+$)','\\1',domain.table$name)
  return(domain.table)
}

#### set WD

args=commandArgs(trailingOnly = T)
#args[1]='/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/ABC_search/AedAeg'
setwd(args[1])

### Import starting data
key=data.table(family=c('ABCA','ABCB','ABCC','ABCD','ABCE','ABCF','ABCG','ABCH','ABC_Unsorted'),
               mins=c(1,1,2,1,2,2,1,1,1),maxes=c(2,2,2,1,2,2,1,1,2))
prelim.fasta=seqinr::read.fasta('Preliminary_ABCs.faa',set.attributes = F,as.string = T,forceDNAtolower = F)


### Parse IPscan
pfam.scan=fread('./HMM_PF00005_output_clean.tsv',select=c(1,3,20,21),sep=' ') %>% rename(name=V1,len=V3,start=V20,end=V21)
pfam.scan$fam=gsub('^.+__(.+)__.+$','\\1',pfam.scan$name)
pfam.table=pfam.scan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% data.table()


#### sort into categories
short.list=list()
good.list=list()
long.list=list()
for(i in unique(pfam.table$fam)){
  sub=pfam.table[fam==i]
  #dom.max=key[family==i]$maxes
  #dom.min=key[family==i]$min
  #short.list[[i]]=sub[domains<dom.min]
  #good.list[[i]]=sub[domains>=dom.min]
  #long.list[[i]]=sub[domains>dom.max]
  short.list[[i]]=sub[len<250]
  good.list[[i]]=sub[len>=250]
  #long.list[[i]]=sub[len>2500]
  
}

### Make tables of each category (short, good, long)
too.short=rbindlist(short.list) %>% domain.annot() %>% data.table()
just.right=rbindlist(good.list) %>% domain.annot() %>% data.table()
#too.long=rbindlist(long.list) %>% domain.annot() %>% data.table()

fwrite(too.short,'./Filtered_short_ABC.tsv')
fwrite(just.right,'./Final_ABC_table.tsv')
#fwrite(too.long,'./Filtered_long_ABC.tsv')


### Make fasta of each category (short, good, long)
too.short.fasta=prelim.fasta[names(prelim.fasta) %in% too.short$name]
just.right.fasta=prelim.fasta[names(prelim.fasta) %in% just.right$name]
#too.long.fasta=prelim.fasta[names(prelim.fasta) %in% too.long$name]

write.fasta(file.out='./Filtered_out_short_ABCS.faa',sequences=too.short.fasta,names=names(too.short.fasta),nbchar=10000)
write.fasta(file.out='./Final_ABCs.faa',sequences=just.right.fasta,names=names(just.right.fasta),nbchar=10000)
#write.fasta(file.out='./Filtered_out_long_ABCS.faa',sequences=too.long.fasta,names=names(too.long.fasta),nbchar=10000)

