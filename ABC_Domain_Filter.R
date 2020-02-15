shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

args = commandArgs(trailingOnly=TRUE)

setwd('/data2/shane/Transporter_ID/ABC_id')

ipscan=fread('./Filter/IPSCAN.tsv',sep='\t',fill=T)[V5=='PF00005'] %>% select(V1,V3,V7,V8) %>% rename(name=V1,len=V3,start=V7,end=V8)
ipscan$fam=gsub('^.+__(.+)__.+$','\\1',ipscan$name)
iptable=ipscan %>% group_by(name,fam) %>% summarize(domains=length(name)) %>% data.table()
blast=fread('./ABC_search/NezVir/total_ABC_recip_blast.tsv') %>% rename(query=V1,subject=V2,evalue=V4,qlen=V5,sstart=V6,send=V7) %>% select(query,subject,evalue,qlen,sstart,send)
blast=blast[!duplicated(query)]


key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted_'),domains=c(2,2,1,2,1,1,1,1,1,1))


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
too.short=rbindlist(short.list)
good=rbindlist(good.list)
too.long=rbindlist(long.list)


domain.annot=function(domain.table,sp){
  test.table=ipscan[name %in% domain.table$name]
  test.table$species=gsub('(^.+)__.+__.+$','\\1',test.table$name)
  test.table$query=gsub('^.+__.+__(.+$)','\\1',test.table$name)
  #test.table[species==sp] 
  blast=fread(paste0('./ABC_search/',sp,'/total_ABC_recip_blast.tsv')) %>% rename(query=V1,subject=V2,evalue=V4,qlen=V5,sstart=V6,send=V7) %>% select(query,subject,evalue,qlen,sstart,send)
  blast=blast[!duplicated(query)] 
  
  final=merge(blast,test.table,by='query') %>% arrange(species,fam,subject) %>% select(-start, -end) %>% unique.data.frame() %>% data.table()
  return(final)
}

nez=domain.annot(too.short,'NezVir')



if(nrow(abc.table)>0){
  colnames(abc.table)=c('code','name')
  cat(format_csv(abc.table))
}
