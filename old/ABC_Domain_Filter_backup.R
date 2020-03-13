#!/usr/bin/env R

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))

#### set up directories and args
setwd('/data2/shane/Transporter_ID/ABC_id')
dir.create('./Filter/preliminary_good',showWarnings=F)
dir.create('./Filter/Short_transporters',showWarnings=F)
dir.create('./Filter/Short_transporters/proteomes',showWarnings=F)
dir.create('./Filter/Short_transporters/dicts',showWarnings=F)
dir.create('./Filter/Long_transporters',showWarnings=F)
dir.create('./Filter/Long_transporters/proteomes',showWarnings=F)
dir.create('./Filter/Long_transporters/dicts',showWarnings=F)
dir.create('./Filter/Counts',showWarnings=F)
dir.create('./Filter/Unsorted_clean',showWarnings=F)

### Import arguments
args = commandArgs(trailingOnly=TRUE)
args[1]=.2
thresh=as.numeric(args[1])

### Import key and metadata
key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted_'),domains=c(2,2,1,2,1,2,2,1,1,1))
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

dict.write=function(x,loc){
	l=x %>% arrange(species) %>% group_split(species)
	sps=sort(x$species) %>% unique()
	for(i in 1:length(l)){
		fwrite(l[[i]],file=paste0('./Filter/',loc,'_transporters/dicts/',sps[i],'.csv'))
	  nam=l[[i]]$name
	  sub.fa=total.fasta[names(total.fasta) %in% nam]
	  write.fasta(sub.fa,names=names(sub.fa),file.out = paste0('./Filter/',loc,'_transporters/proteomes/',sps[i],'.faa'))
	}
}


fasta.write=function(x,loc){
  l=x %>% arrange(species) %>% group_split(species)
  sps=sort(x$species) %>% unique()
  for(i in 1:length(l)){
    fwrite(l[[i]],file=paste0('./Filter/',loc,'_transporters/dicts/',sps[i],'.csv'))
  }
}


### Parse IPscan
#ipscan=fread('/data2/shane/Transporter_ID/ABC_id/Filter/old_IP/IP_combined.tsv',sep='\t',fill=T)[V5=='PF00005'] %>% select(V1,V3,V7,V8) %>% rename(name=V1,len=V3,start=V7,end=V8)
ipscan=fread('./Filter/IPSCAN.tsv',sep='\t',fill=T)[V5=='PF00005'] %>% select(V1,V3,V7,V8) %>% rename(name=V1,len=V3,start=V7,end=V8)
ipscan$fam=gsub('^.+__(.+)__.+$','\\1',ipscan$name)
iptable=ipscan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% data.table()
total.fasta=read.fasta('./Filter/ABC_preliminary_total.faa',set.attributes = F,as.string = T,forceDNAtolower = F)

short.list=list()
good.list=list()
long.list=list()
for(i in unique(iptable$fam)){
  sub=iptable[fam==i]
  mins=key[family==i]$domains
  short.list[[i]]=sub[domains<mins]
  good.list[[i]]=sub[domains==mins]
  long.list[[i]]=sub[domains>mins]
}
too.short=rbindlist(short.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()
good=rbindlist(good.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()
too.long=rbindlist(long.list) %>% filter(!grepl('Unsorted',name)) %>% domain.annot() %>% data.table()

#### Deal with unsorted proteins
unsorted.good=iptable[fam=='ABC_Unsorted_' & len >100 & domains<3]
unsorted.fa=total.fasta[(names(total.fasta) %in% unsorted.good$name) | (grepl('DroMel',names(total.fasta)))]
write.fasta(sequences = unsorted.fa,names = gsub(':','_',names(unsorted.fa),fixed=T),file.out = './Filter/Unsorted_clean/Unsorted_good.faa')


###### ADD in model species to good

#### WRITE dictionaries
dict.write(too.short,'Short')
dict.write(too.long,'Long')
dict.write(good,'Full')

#### produce counts files
good.sum=good %>% count.fams() %>% spread(key=species,value=count) %>%
     data.table()
good.sum[is.na(good.sum)]=0
good.sum=data.table(fam=good.sum$fam,apply(select(good.sum,2:length(colnames(good.sum))),2,as.numeric)) 
fwrite(good.sum,'./Filter/Counts/Full_transporter_counts.csv') 


#apply(select(good.sum,2:length(colnames(good.sum))),2,min)

good.count.sp=good %>% count.fams() %>% group_by(species) %>% summarize(good.count=sum(count)) %>% data.table()
good.count.fam=good %>% count.fams() %>% data.table()

short.count.sp=too.short %>% count.fams() %>% group_by(species) %>% summarize(short.count=sum(count)) %>% data.table()
short.count.fam=too.short %>% count.fams() %>% data.table()

long.count.sp=too.long %>% count.fams() %>% group_by(species) %>% summarize(long.count=sum(count)) %>% data.table()
long.count.fam=too.long %>% count.fams() %>% data.table()

m=merge(good.count.sp,short.count.sp,by='species',all=T) %>% merge(long.count.sp,by='species',all=T) %>% 
  mutate(total.count=good.count+short.count+long.count,per_short=short.count/total.count) %>% data.table()
m[is.na(m)]=0
good.species=m[per_short<thresh]$species

quality_cutoff=metadata[abbreviation %in% good.species]

final.good.species=quality_cutoff$taxid_code
writeLines(final.good.species,'./Filter/Quality_cutoff_species.txt')
fwrite(quality_cutoff,'./Filter/Quality_cutoff_species.tsv',sep='\t')

fwrite(m,'./Filter/Counts/short_long_summary')


fwrite(short.count.sp,'./Filter/Counts/short_summary_species.csv')
fwrite(short.count.fam,'./Filter/Counts/short_summary_family.csv')


fwrite(long.count.sp,'./Filter/Counts/long_summary_species.csv')
fwrite(long.count.fam,'./Filter/Counts/long_summary_family.csv')


