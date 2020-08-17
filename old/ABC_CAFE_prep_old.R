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

#setwd('/mnt/disk/shane/Transporter_ID/ABC_id')
dir.create('./CAFE/CAFE_tables',showWarnings = F)

### Import data
abc.counts=fread('./Final_outputs/combined_files/Full_counts_wide.tsv')
used.species=colnames(abc.counts)[2:length(colnames(abc.counts))]

#metadata=fread('./GENERAL_REFERENCE/keys/Arthropod_species_metadata.tsv',header=T) %>% 
#  select(Species_name,abbreviation,taxid_code) %>%  data.table()


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
fwrite(abc.counts,'./CAFE/ABC_COUNTS_CAFE_FULL.tsv',sep='\t')


### create list of tree base names
iter=list.files('./CAFE/clean_raxml_trees')[grepl('RAxML_bipartitions.',list.files('./CAFE/clean_raxml_trees'))] %>%
  str_remove('RAxML_bipartitions.') %>% str_remove('_species.nwk')



###### RUN LOOP
for (i in iter){
  
  ### rename tree
  tr=read.tree(paste0('./CAFE/clean_raxml_trees/RAxML_bipartitions.',i,'_species.nwk'))
  
  #### this chunk renmaes taxid codes to 6 letter abbreviations. It has been deprecated with the latest species tree phylogeny script.
  #for(j in metadata$taxid_code){
  #  if(grepl(j,rax.tree)){
  #    rax.tree=gsub(j,metadata[taxid_code==j]$abbreviation,rax.tree)
  #  }
  #}
  #writeLines(rax.tree,paste0('./CAFE/clean_raxml_trees/raxml_tree_named_',i,'.tre'))
  
  ## read in named raxmL tree
  #tr=read.tree(paste0('./CAFE/clean_raxml_trees/raxml_tree_named_',i,'.tre'))
  
  ### Remove tips which are not in final dataset
  drops=tr$tip.label[!(tr$tip.label %in% used.species)]
  tr=drop.tip(tr,drops)
  #tr=drop.tip(tr,tr$tip.label[grepl('[0-9]+_0',tr$tip.label)])
  #if(i=='Diptera'){tr=drop.tip(tr,'BomMor')}
  if(grepl('ArachInsect',i)){tr=root(tr,outgroup=c('LepDel','TetUrt','SteMim','ParTep','CenScu','TroMer','VarDes','GalOcc','IxoSca'))}
  
  ### Import calibrations for various nodes
  nodes <- c(); maxes=c()
  maxes=c()
  mins=c()
  if(("AcyPis" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","DroMel")));maxes=c(maxes,401);mins=c(mins,345)}
  if(("ApiMel" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ApiMel","DroMel")));maxes=c(maxes,372);mins=c(mins,317)} ## has fossil
  #if(("AedAeg" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AedAeg","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
  if(("AnoGam" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AnoGam","DroMel")));maxes=c(maxes,206);mins=c(mins,107)} ## has fossil
  if(("NilLug" %in% tr$tip.label) & ("AcyPis" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("AcyPis","NilLug")));maxes=c(maxes,346);mins=c(mins,232)}
  #if(("TetUrt" %in% tr$tip.label) & ("DroMel" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("TetUrt","DroMel")));maxes=c(maxes,579);mins=c(mins,539)}
  if(("PluXyl" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("PluXyl","BomMor")));maxes=c(maxes,178);mins=c(mins,116)} ## has fossil
  if(("DroMel" %in% tr$tip.label) & ("BomMor" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("DroMel","BomMor")));maxes=c(maxes,328);mins=c(mins,244)} ## has fossil  ## create ultrametric tree
  if(("IxoSca" %in% tr$tip.label) & ("VarDes" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("IxoSca","VarDes")));maxes=c(maxes,384);mins=c(mins,282)} 
  #if(("ParTep" %in% tr$tip.label) & ("CenScu" %in% tr$tip.label)){nodes=c(nodes,getMRCA(tr, tip = c("ParTep","CenScu")));maxes=c(maxes,419);mins=c(mins,354)} 
  
   
  ## arachind dates from https://link.springer.com/article/10.1007/s10493-008-9203-5/tables/4
  
  
  ### Credit to Alex SL for format https://phylobotanist.blogspot.com/2019/
  mycalibration <- makeChronosCalib(tr, node=c(nodes), age.min=mins,age.max=maxes)
  mytimetree <- chronos(tr, lambda = 1, model = "correlated", calibration = mycalibration)
  num=mytimetree$node.label %>% as.numeric()
  write.tree(mytimetree, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))
  
  ### write ultrametric tree without nodelabels for CAFE processing. CAFE python script will not work without this
  mytimetree.nolab=mytimetree
  mytimetree.nolab$node.label=NULL
  write.tree(mytimetree.nolab, file=paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric_nolab.tre'))
  
  #### Prepare tree to plot
  plot.tree=read.tree(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))
  cols=c()
  for(j in num){
    if(is.na(j)){cols=c(cols,'black')
    }else if(j>80){cols=c(cols,'green4')
    }else if(j>50){cols=c(cols,'blue4')
    }else{cols=c(cols,'red')}
  }
  ma=max(mytimetree$edge.length)
  xma=ma+75
  ma.r=seq(0,round(ma,-2),by=100)
  diff=ma-round(ma,-2)
  
  ### make plot 
  gp=ggtree(plot.tree)#, mrsd = "2010-01-01")
  gp=gp+geom_tiplab(size=6)
  gp=gp+geom_nodepoint(size=5,col=cols)
  gp=gp+xlab('Millions of years ago (Mya)')
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  
  ggsave(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.pdf'),plot=gp,width=14,height=10)
  
  
  ## create lambda file
  writeLines(lambda.convert(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre'))),paste0("./CAFE/clean_raxml_trees/",i,'_tree_lambda.txt')) 
  
  sp=str_extract_all(readLines(paste0("./CAFE/clean_raxml_trees/",i,'_tree_ultrametric.tre')),pattern = "[A-z]+",simplify = T)  %>% as.character()
  
  sp.counts=abc.counts %>% select(c('Desc','Family ID',all_of(sp)))
  fwrite(sp.counts,paste0('./CAFE/CAFE_tables/',i,'_ABC_CAFE_table.tsv'),sep='\t')
}





