#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))
shhh(library(ggplot2))

#### FILTER OUT QUALITY ON LINE 61

setwd('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline')
dir.create('./CAFE/CAFE_tables',showWarnings = F)

### Import data
abc.counts=fread('./Final_outputs/combined_files/Full_counts_wide.tsv')
used.species=fread('./Final_outputs/combined_files/Full_counts_long.tsv')$abbreviation
#colnames(abc.counts)[2:length(colnames(abc.counts))]

full.metadata=fread('./GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv',header=T) %>% 
 select(Species_name,abbreviation,taxid_code) %>%  data.table()


### Import functions
lambda.convert=function(x){
  a=gsub("[A-z]+:[0-9]+.[0-9]+",1,x)
  b=gsub("[0-9]{2,}:[0-9]+.[0-9]+",1,a)
  c=gsub(":[0-9]+.[0-9]+",1,b)
  d=gsub('[A-z]+:[0-9]+',1,c)
  return(d)
}

### format table for CAFE
colnames(abc.counts)[1]='Family ID'
abc.counts$Desc='(null)'
abc.counts=select(abc.counts,Desc,'Family ID',everything())
abc.counts=abc.counts[`Family ID`!='ABC_Unsorted']
fwrite(abc.counts,'./CAFE/ABC_COUNTS_CAFE_FULL.tsv',sep='\t')

###calculate xlab functions
xma.calc=function(tree){
  table=as_tibble(tree) %>% data.table()
  sum1=table[node %in% tidytree::ancestor(tree,1)]$branch.length %>% sum(na.rm=T)
  sum2=table[node==1]$branch.length
  return(sum1+sum2)
}
  
### create list of tree base names
iter=list.files('./CAFE/clean_raxml_trees')[grepl('Raw.+nwk',list.files('./CAFE/clean_raxml_trees'))] %>%
  str_remove('_species_tree.nwk') %>% str_remove('Raw_')

###### RUN LOOP
for (i in iter){
  
  tr1=read.tree(paste0('./CAFE/clean_raxml_trees/Raw_',i,'_species_tree.nwk'))
  drops=tr1$tip.label[!(tr1$tip.label %in% used.species)]
  if(grepl('Arachnid',i)){tr1=root(tr1,outgroup=c('SarSca','LepDel','TetUrt','SteMim','ParTep','CenScu','TroMer','VarAac','VarDes','GalOcc','IxoSca'))}
  if(grepl('Coleoptera',i)){tr1=root(tr1,outgroup=c('ApiMel'))}
  if(grepl('Diptera',i)){
  tr1=root(tr1,outgroup=c('BomMor'))
  drops=c(drops,'BomMor')
  }
  if(grepl('Hemimetabola',i)){tr1=root(tr1,outgroup=c('DroMel','ApiMel','AedAeg','BomMor'))}
  if(grepl('Lepidoptera',i)){
  tr1=root(tr1,outgroup=c("ApiMel"))
  #drops=c(drops,'MyzPer')
  }
  if(grepl('Hymenoptera',i)){
    tr1=root(tr1,outgroup=c('AthRos','CepCin','NeoLec'))
    #drops=c(drops,'MyzPer')
  }
  
  
  
  tr=drop.tip(tr1,drops)
  
  
  write.tree(tr,paste0('./CAFE/clean_raxml_trees/Clean_',i,'_species_tree.nwk'))
  
  nodes <- c(); maxes=c()
  maxes=c()
  mins=c()
  #if(("CaeEle" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("CaeEle","DroMel")));maxes=c(maxes,1000);mins=c(800)}
  if(("AcyPis" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,401);mins=c(mins,345)}
  if(("ApiMel" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ApiMel","DroMel")));maxes=c(maxes,372);mins=c(mins,317)} ## has fossil
  if(("AedAeg" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AedAeg","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
  if(("NilLug" %in% tr$tip.label) & ("AcyPis" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","NilLug")));maxes=c(maxes,346);mins=c(mins,232)}
  if(("TetUrt" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("TetUrt","DroMel")));maxes=c(maxes,579);mins=c(mins,539)}
  if(("PluXyl" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("PluXyl","BomMor")));maxes=c(maxes,178);mins=c(mins,116)} ## has fossil
  if(("DroMel" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("DroMel","BomMor")));maxes=c(maxes,328);mins=c(mins,244)} ## has fossil
  if(("IxoSca" %in% tr$tip.label) & ("VarDes" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("IxoSca","VarDes")));maxes=c(maxes,384);mins=c(mins,282)} 
  if(("ParTep" %in% tr$tip.label) & ("CenScu" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ParTep","CenScu")));maxes=c(maxes,419);mins=c(mins,354)} 
  if(("BomTer" %in% tr$tip.label) & ("ApiMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("BomTer","ApiMel")));maxes=c(maxes,40);mins=c(mins,20)} 
  if(("AcrEch" %in% tr$tip.label) & ("ApiMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcrEch","ApiMel")));maxes=c(maxes,110);mins=c(mins,90)} 
  if(("OruAbi" %in% tr$tip.label) & ("ApiMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("OruAbi","ApiMel")));maxes=c(maxes,180);mins=c(mins,140)} 
  if(("NeoLec" %in% tr$tip.label) & ("ApiMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("NeoLec","ApiMel")));maxes=c(maxes,260);mins=c(mins,240)} 
  
  
  
  ## create ultrametric tree
  ### Credit to Alex SL for format https://phylobotanist.blogspot.com/2019/
  mycalibration <- makeChronosCalib(tr, node=c(nodes), age.min=mins,age.max=maxes)
  mytimetree <- chronos(tr, lambda = 1, model = "correlated", calibration = mycalibration)
  num=mytimetree$node.label %>% as.numeric()
  mytimetree$node.label=NULL
  #mytimetree$edge.length=round(mytimetree$edge.length)
  write.tree(mytimetree, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.nwk'))
  
  
  #### Plot tree
  plot.tree=read.tree(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.nwk'))
  cols=c()
  for(j in num){
    if(is.na(j)){cols=c(cols,'black')
    }else if(j>80){cols=c(cols,'green4')
    }else if(j>50){cols=c(cols,'blue4')
    }else{cols=c(cols,'red')}
  }
  #plot.tree$node
  ma=xma.calc(plot.tree)
  #ma=max(mytimetree$edge.length)
  xma=ma+100
  ma.r=seq(0,round(ma,-2),by=100)
  
  diff=ma-round(ma,-2)
  
  #pdf(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.pdf'),width=14,height=7)
  gp=ggtree(plot.tree)#, mrsd = "2010-01-01")
  gp=gp+geom_tiplab(size=6)
  gp=gp+geom_nodepoint(size=5,col=cols)
  gp=gp+xlab('Millions of years ago (Mya)')
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(-50,xma))
  #print(gp)
  #dev.off()
  
  ggsave(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.pdf'),plot=gp,width=14,height=10)
  
  #l.nwke.ch=chronopl(read.nwke(paste0(paste0(H,'CAFE/trees/raxml_tree_named_',i,'.nwk')), lambda=0.1)
  #l.nwke.ch$edge.length=l.nwke.ch$edge.length*1000
  #l.nwke.ch$node.label=NULL
  #write.nwke(l.nwke.ch, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.nwk'))
  
  ## create lambda file
  writeLines(lambda.convert(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.nwk'))),paste0("./CAFE/clean_raxml_trees/",i,'_tree_lambda.txt')) 
  
  
  sp=str_extract_all(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.nwk')),pattern = "[A-z]+",simplify = T)
  sp=as.character(sp)
  
  sp.counts=abc.counts %>% select(c('Desc','Family ID',all_of(sp)))
  fwrite(sp.counts,paste0('./CAFE/CAFE_tables/',i,'_ABC_CAFE_table.tsv'),sep='\t')
  #orthodb.counts=orthodb.final %>% select(c('Desc','Family ID',sp))
  #fwrite(orthodb.counts,paste0('./CAFE/CAFE_tables/',i,'_OrthoDB_CAFE_table.tsv'),sep='\t')
}




