#!/usr/bin/env Rscript

shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))
shhh(library(seqinr))
shhh(library(argparser))


#### Import arguments
p=arg_parser('Combine ABC table')
p = add_argument(p, "--minlen", help="minimum length",default=250)
p = add_argument(p, "--motif",help='must have this domain in sequence',default=NA)
p = add_argument(p, "--domain_filter",help='use domain filter key',default=NA)
p = add_argument(p, "--indir",help='use domain filter key',default='/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/ABC_search/DanPle')

argv=parse_args(p)

minlen=argv$minlen
motif=argv$motif
domain_filter=argv$domain_filter
indir=argv$indir


if(!is.na(motif)){
  if(motif=="NA"){
    motif=NA
  }
}

if(!is.na(domain_filter)){
  if(domain_filter=="NA"){
    domain_filter=NA
  }
}

#### set WD
setwd(indir)


### Functions
domain.annot=function(domain.table){
  domain.table$abbreviation=gsub('(^.+)__.+__.+$','\\1',domain.table$name)
  domain.table$code=gsub('^.+__.+__(.+$)','\\1',domain.table$name)
  return(domain.table)
}


### Import starting data
key=data.table(family=c('ABCA','ABCB','ABCC','ABCD','ABCE','ABCF','ABCG','ABCH','ABC_Unsorted'),
               mins=c(1,1,2,1,2,2,1,1,1),maxes=c(2,2,2,1,2,2,1,1,2))
prelim.fasta=seqinr::read.fasta('Preliminary_ABCs.faa',set.attributes = F,as.string = T,forceDNAtolower = F)

### Parse IPscan
pfam.scan=fread('./HMM_PF00005_output_clean.tsv',select=c(1,3,20,21),sep=' ') %>% rename(name=V1,len=V3,start=V20,end=V21)
pfam.scan$fam=gsub('^.+__(.+)__.+$','\\1',pfam.scan$name)
pfam.table=pfam.scan %>% group_by(name,fam,len) %>% summarize(domains=length(name)) %>% data.table()


## Add in motifs
motifs=fread('./Mem_motif.tsv',header=F,sep='\t')
motif.clean=dcast(motifs,V1~V2,value.var='V2',fun.aggregate = function(x) ifelse(length(x) > 0,'Present','Absent'))
colnames(motif.clean)=c('name','WalkerB','WalkerA','ABC_sig')


total.table=merge(pfam.table,motif.clean,by='name',all=T)
total.table[is.na(total.table)]='Absent'


#### sort into categories
filt.list=list()
good.list=list()
for(i in total.table$name){
  sub=total.table[name==i]
  dom.min=key[family==sub$fam]$min
  dom.max=key[family==sub$fam]$max
  
  if(sub$len<minlen){
    filt.list[[i]]=sub
  }else if (!is.na(domain_filter)){
    if(!(sub$domains>=dom.min & sub$domains<=dom.max)){
      filt.list[[i]]=sub
    }
  }else if(motif %in% colnames(sub)){
    if(sub[[motif]]=='Absent'){
      filt.list[[i]]=sub
    }
  }else{
    good.list[[i]]=sub
  }
}

### Make tables of each category (short, good, long)
filt=rbindlist(filt.list) %>% domain.annot() %>% data.table()
good=rbindlist(good.list) %>% domain.annot() %>% data.table()
total=rbind(filt,good,fill=T)
#too.long=rbindlist(long.list) %>% domain.annot() %>% data.table()

fwrite(filt,'./Filtered_ABC.tsv')
fwrite(good,'./Final_ABC_table.tsv')
fwrite(total,'./Total_annotated_ABC.tsv')


### Make fasta of each category (short, good, long)
filt.fasta=prelim.fasta[names(prelim.fasta) %in% filt$name]
good.fasta=prelim.fasta[names(prelim.fasta) %in% good$name]
total.fasta=prelim.fasta[names(prelim.fasta) %in% total$name]


write.fasta(file.out='./Filtered_ABCS.faa',sequences=filt.fasta,names=names(filt.fasta),nbchar=10000)
write.fasta(file.out='./Final_ABCs.faa',sequences=good.fasta,names=names(good.fasta),nbchar=10000)
write.fasta(file.out='./Total_annotated_ABCs.faa',sequences=total.fasta,names=names(total.fasta),nbchar=10000)

