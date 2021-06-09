#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(tidyverse))
shhh(library(data.table))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(treeio))
setwd('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/')
dir.create('./CAFE/CAFE_figures')

used.species=fread('./Final_outputs/combined_files/Full_counts_long.tsv')$abbreviation


iter=list.files('./CAFE/CAFE_tables/') %>% gsub('_ABC_CAFE_table.tsv','',.)
full.metadata=fread('./Final_outputs/combined_files/Full_counts_long.tsv')
#group='Coleoptera_taxid_codes'
group='Lepidoptera'
family='ABCH'
node.annot=''
label.annot=''
nodes=F

short.spec=function(x){
  gen=substr(unlist(str_split(x,'_'))[1],1,1)
  specific=unlist(str_split(x,'_'))[2]
  final=paste0(gen,'. ',specific)
}


xma.calc=function(tree){
  table=as_tibble(tree) %>% data.table()
  sum1=table[node %in% tidytree::ancestor(tree,1)]$branch.length %>% sum(na.rm=T)
  sum2=table[node==1]$branch.length
  return(sum1+sum2)
}

### tree.fig function                 
tree.fig=function(group,family,node.annot='',label.annot='',nodes=F){
  
  ##import ultrametric tree
  ultra.tree=read.tree(paste0('./CAFE/clean_raxml_trees/',group,'_tree_ultrametric.nwk'))
  tbl=as_tibble(ultra.tree) %>% data.table()
  
  #### set colors
  node.tree=read.tree(paste0('./CAFE/clean_raxml_trees/Clean_',group,'_species_tree.nwk'))
  drops=node.tree$tip.label[!(node.tree$tip.label %in% used.species)]
  node.tree=drop.tip(node.tree,drops)

  node.scores=as.numeric(node.tree$node.label)
  cols=c()
  for(j in node.scores){
    if(is.na(j)){cols=c(cols,'grey50')
    }else if(j>80){cols=c(cols,'green4')
    }else if(j>50){cols=c(cols,'gold4')
    }else{cols=c(cols,'red')}
  }
  
  if(ultra.tree$Nnode>length(cols)){
  	cols=c('grey50',cols)
  }else if(ultra.tree$Nnode<length(cols)){
  	cols=cols[2:length(cols)]
  }
  
  
  
  ### Import node labels
  lab.text=readLines(paste0('./CAFE/outputs/',group,'/','Base_asr.tre'))[3] %>% 
    gsub("^.+ = ",'',.) %>% gsub('>_[0-9|\\:|\\.]+','>',.) %>%
    gsub(';','',.) 
  lab.text=gsub('\\*_[0-9]+','',lab.text)

  lab.tree=read.tree(text=paste0('(',lab.text,')',';'))
    
  ## import counts
  #count.table=fread('/mnt/disk/shane/Transporter_ID/SLC_id_pipeline/CAFE/outputs/ArachInsect_species/Base_count.tab')[`FamilyID`==family]
  count.table=fread(paste0('./CAFE/outputs/',group,'/','Base_count.tab'))[`FamilyID`==family]
  #count.table=fread(paste0('./CAFE/outputs/',group,'/','Base_count.tab')) %>% rename(family=`FamilyID`)
  if(nrow(count.table)==0){
    break
  }
  count.reduce=count.table %>% select(-matches('[A-z]'))
  count.term=count.table %>% select(matches('[A-z]')) %>% select(-`FamilyID`)
  colnames(count.term)=gsub('<[0-9]+>','',colnames(count.term))
  
  ## add node labels to ultrametric
  sorted=c()
  for(i in lab.tree$node.label){
    temp=count.reduce[[i]]
    sorted=c(sorted,temp)
  }
  final=c('',as.character(sorted))
  final2=final[-1]
  
  
  ultra.tree$node.label=final2
  if(nodes==T){ultra.tree$node.label=lab.tree$node.label[-1]} ## whether to create tree with node labels instead of CAFE counts
  
  ## add tip labels to ultrametric with numbers
  for(i in ultra.tree$tip.label){
    temp=count.term[[i]]
    species.name=short.spec(full.metadata[abbreviation==i]$Species_name)
    ultra.tree$tip.label[which(ultra.tree$tip.label==i)]=paste0(species.name,' (',temp,')')
  }
  
  ##Set scaling factors
  ma=xma.calc(ultra.tree)
  xma=ma*1.5
  ma.r=seq(0,round(ma,-2),by=100)
  diff=ma-round(ma,-2)
  
  
  
  #### Make plot
  
  gp=ggtree(ultra.tree,size=2)
  
  if(grepl('Diptera',group)){
  gp=gp+geom_tiplab(size=6,fontface='bold')#,aes(label=paste0('bold(', label, ')')), parse=TRUE)
  gp=gp+geom_nodepoint(size=10,color=cols)
  gp=gp+geom_nodelab(hjust=.6,size=4,fontface='bold',color='white')
  
  }else{
 
  gp=gp+geom_tiplab(size=12,fontface='bold')#,aes(label=paste0('bold(', label, ')')), parse=TRUE)
  gp=gp+geom_nodepoint(size=21,color=cols)
  gp=gp+geom_nodelab(hjust=.6,size=12,fontface='bold',color='white')
  }
  if(is.list(node.annot)){
    names(node.annot)=label.annot
    plot.annot=vector("list",length=length(label.annot))
    names(plot.annot)=label.annot
    for(i in names(plot.annot)){
      plot.annot[[i]][1]=tbl[grepl(node.annot[[i]][1],label)]$node %>% as.numeric()
      plot.annot[[i]][2]=tbl[grepl(node.annot[[i]][2],label)]$node %>% as.numeric()
      gp=gp+geom_strip(taxa1=plot.annot[[i]][1],taxa2=plot.annot[[i]][2],offset.text=4,fontsize=10,
                       barsize=5,color='black',label=i,offset=ma/6.5)
    }
  }
  #gp=gp+geom_hilight(node=list(node1,55),fill='darkgreen',alpha=.3)
  gp=gp+labs(x='Millions of Years Ago')
  #gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=30,face='bold',color = 'black'),axis.line.x=element_line(size=3),
              axis.title.x=element_text(size=40,face='bold'),
              plot.margin=unit(c(t=0,r=2,b=0,l=0),"cm"))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  #print(gp)
  #ggsave('test.pdf',gp,width=16,height=9)
}

a=tree.fig(group='Coleoptera',family='ABCA')
fams=colnames(full.metadata)[grepl('ABC',colnames(full.metadata))]
fams=fams[fams!='SLC_14' & fams!='ABC_Unsorted' & fams!='ABC_total']
  

#iter=iter[2]
for (i in iter){ 
  counts=fread(paste0('./CAFE/CAFE_tables/',i,'_ABC_CAFE_table.tsv'))
  for(j in fams){ 
    res=try(tree.fig(group = i,family = j),silent=T)
    if(class(res)!='try-error'){
      temp=tree.fig(group = i,family = j)
      ggsave(plot=temp,filename=paste0('./CAFE/CAFE_figures/',i,'_',j,'.pdf'),device='pdf',width=16,height=11)
    }
  } 
} 
