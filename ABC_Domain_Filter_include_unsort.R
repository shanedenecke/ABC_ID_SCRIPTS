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
dir.create('./Filter/Final_transporters')
dir.create('./Filter/Final_transporters/dicts')
dir.create('./Filter/Final_transporters/proteomes')
dir.create('./Filter/Total_combined')


key=data.table(family=c(gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),'ABC_Unsorted_'),domains=c(2,2,1,2,1,2,2,1,1,1))
metadata=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv',header=T) %>% 
  select(Species_name,abbreviation,taxid_code)
good.species=readLines('./Filter/preliminary/Quality_threshold_species.txt')
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

##############################################
########## read and process unsorted tree
#unsorted.tree=read.tree('./Filter/Unsorted_clean/RAxML_bipartitions.Unsorted.nwk')
#collapsed.tree=di2multi4node(unsorted.tree,10)
### mark drosophila as red
#cols=rep('black',length=length(unsorted.tree$tip.label))
#cols[grepl('DroMel',unsorted.tree$tip.label)]='red'
#cols[grepl('TriCas',unsorted.tree$tip.label)]='blue'
#gp=ggtree(unsorted.tree,size=2,layout='rectangular')
#gp=gp+geom_tiplab(size=4,fontface='bold',color=cols)
#gp=gp+geom_nodepoint(size=4,col='black')
#gp=gp+geom_nodelab(hjust=1,vjust=.3,size=1.5,fontface='bold',col='white')
#gp=gp+lims(x=c(0,20))
#gp=gp+theme(title = element_text(size=12))
#print(gp)
#ggsave(plot=gp,filename ='./Filter/Unsorted_clean/Unosrted_ABC_phylogeny.pdf',device='pdf',height=40,width=30,limitsize = F)
#ggsave(plot=gp,filename ='~/Dropbox/Unosrted_ABC_phylogeny.pdf',device='pdf',height=40,width=30,limitsize = F)


#tbl=as_tibble(unsorted.tree)
#dros.tips=unsorted.tree$tip.label[grepl('DroMel',unsorted.tree$tip.label)]
#fams=dros.tips %>% gsub('DroMel__(ABC[A-Z]+)_.+$','\\1',.) %>% unique()


#### manual
#a.node=MRCA(tbl,dros.tips[grepl('ABCA',dros.tips)])$node
#f.node=MRCA(tbl,'DroMel__ABCF__FBgn0030672_CG9281','CimLec__ABC_Unsorted___XP_0142575901')$node
#g.node=MRCA(tbl,'DroMel__ABCG__FBgn0052091_CG32091','CimLec__ABC_Unsorted___XP_0142577721')$node
#b.node=MRCA(tbl,'DroMel__ABCBF__FBgn0004513_Mdr65','HyaAzt__ABC_Unsorted___LOC108683648')$node
#d.node=MRCA(tbl,'DroObs__ABC_Unsorted___LOC111067316','MenMol__ABC_Unsorted___1155016_002EAC')$node

#l=list()
#for(i in fams){
#  sub=dros.tips[grepl(i,dros.tips)]
#  fam.node=parent(tbl,MRCA(tbl,sub)$node)$node
#  fam.tbl=offspring(tbl,fam.node) %>% filter(!grepl('DroMel',label) & grepl('[A-Z]',label))
#  fam.genes=fam.tbl$label
#  ip.sub=iptable[name %in% fam.genes]
#  l[[i]]=ip.sub
#}
#unsorted.sum=rbindlist(l)


############################################################
#NEED TO CREATE IPTABLE WITH UNSORTED 
# ALSO NEED TO RENAME FASTA FILE




good.list=list()
for(i in unique(iptable$fam)){
  sub=iptable[fam==i]
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
