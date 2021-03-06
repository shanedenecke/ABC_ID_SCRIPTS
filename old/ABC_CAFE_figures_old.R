shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(tidyr))
shhh(library(ggplot2))
shhh(library(treeio))

### set directories
#setwd('/mnt/disk/shane/Transporter_ID/ABC_id/')
dir.create('./CAFE/CAFE_figures',showWarnings = F)

#### Write CAFE count tables for each tre What are these used for?
iter=list.files('./CAFE/CAFE_tables/') %>% gsub('_ABC_CAFE_table.tsv','',.)
counts=fread('./CAFE/ABC_COUNTS_CAFE_FULL.tsv') %>% select(-Desc) %>% rename(fam=`Family ID`)


for(i in iter){
  ids=gsub("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ",
           "",readLines(paste0('./CAFE/outputs/',i,'_ABC_cafe_output.cafe'))[5],fixed=T)
  ids2=unlist(strsplit(ids,split=' '))
  sp=gsub('# IDs of nodes:','',readLines(paste0('./CAFE/outputs/',i,'_ABC_cafe_output.cafe'))[4])
  sp2=unlist(str_extract_all(sp,"[A-z]{6,}<[0-9]+>")) %>% str_remove_all('<|>')
  cafe=fread(paste0('./CAFE/outputs/',i,'_ABC_cafe_output.cafe'),skip=11,sep='\t') %>% select(V1,V3,V4) %>% 
    filter(V3<.05) %>% separate(col=V4,into=ids2,sep='\\),\\(') %>%
    data.table()
  fwrite(cafe,paste0('./CAFE/CAFE_figures/',i,'_CAFE.csv'))
}
  

#group='Diptera'
#family='ABCA'
#node.annot = list(c('DroEre','DroBus'),
#                  c('CluMar','AnoCul'))

#label.annot=c('Drosophilid','Mosquito')          
 

### tree.fig function                 
tree.fig=function(group,family,node.annot='',label.annot=''){
  
  ##import ultrametric tree
  base.tree=read.tree(paste0('./CAFE/clean_raxml_trees/',group,'_tree_ultrametric.tre'))
  tbl=as_tibble(base.tree) %>% data.table()
  
  #### set colors
  node.scores=as.numeric(base.tree$node.label)
  cols=c()
  for(j in node.scores){
    if(is.na(j)){cols=c(cols,'grey50')
    }else if(j>90){cols=c(cols,'green4')
    }else if(j>70){cols=c(cols,'gold4')
    }else{cols=c(cols,'red')}
  }
  
  
  ### Import node labels
  lab.text=gsub('# The labeled CAFE tree:\t','',readLines(paste0('./CAFE/outputs/',group,'_ABC_summary.txt_fams.txt'))[1])
  lab.tree=read.tree(text=paste0('(',lab.text,')',';'))
  
  ## import counts
  count.table=fread(paste0('./CAFE/outputs/',group,'_ABC_summary.txt_anc.txt'))[`Family ID`==family]
  count.reduce=count.table %>% select(-matches('[A-z]'))
  count.term=count.table %>% select(matches('[A-z]')) %>% select(-`Family ID`)
  colnames(count.term)=gsub('<[0-9]+>','',colnames(count.term))
  
  ## add node labels to ultrametric
  sorted=c()
  for(i in lab.tree$node.label){
    temp=count.reduce[[i]]
    sorted=c(sorted,temp)
  }
  final=c('',as.character(sorted))
  base.tree$node.label=final
  
  ## add tip labels to ultrametric with numbers
  for(i in base.tree$tip.label){
    temp=count.term[[i]]
    base.tree$tip.label[which(base.tree$tip.label==i)]=paste0(i,' (',temp,')')
  }
  
  ##Set scaling factors
  #ma=max(sapply(1:base.tree$Nnode,function(x) tbl[node %in% ancestor(base.tree,x)]$branch.length %>% sum(na.rm = T)))
  #xma=ma*1.4
  #ma.r=seq(0,round(ma,-2),by=50)
  #diff=ma-round(ma,-2)
  
  ma=max(base.tree$edge.length)
  xma=ma*1.3
  ma.r=seq(0,round(ma,-2),by=100)
  diff=ma-round(ma,-2)
  
  
  
  #### Make plot
  
  gp=ggtree(base.tree,size=2)
  gp=gp+geom_tiplab(size=8,fontface='bold')#,aes(label=paste0('bold(', label, ')')), parse=TRUE)
  gp=gp+geom_nodepoint(size=16,color=cols)
  gp=gp+geom_nodelab(hjust=.75,size=8,fontface='bold',color='white')
  
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
  gp=gp+theme_tree2()
  gp=gp+theme(axis.text.x=element_text(size=20,face='bold',color = 'black'),axis.line.x=element_line(size=3),
              axis.title.x=element_text(size=20))
  gp=gp+scale_x_continuous(breaks=diff+ma.r,labels=as.character(rev(ma.r)),limits=c(0,xma))
  print(gp)
}


for (i in iter){ 
  for(j in counts$fam[!grepl('Unsorted',counts$fam)]){ 
    temp=tree.fig(group = i,family = j)
    ggsave(plot=temp,filename=paste0('./CAFE/CAFE_figures/',i,'_',j,'.pdf'),device='pdf',width=25,height=15)
  }
}



#dros.a=tree.fig(group = 'Diptera',family = 'ABCA',node.annot = list(c('DroEre','DroBus'),c('CluMar','AnoCul')),
#         label.annot=c('Drosophilid','Mosquito'))
#ggsave(plot=dros.a,filename='./CAFE/CAFE_figures/Drosophila_ABCA.pdf',device='pdf',width=25,height=15)


#hemi.h=tree.fig(group = 'Hemipteran',family = 'ABCH',node.annot = list(c('BomMor','DroMel'),c('MyzPer','AphGly')),
#         label.annot=c('Holometabola','Aphid'))
#ggsave(plot=hemi.h,filename='./CAFE/CAFE_figures/Hemi_H.pdf',device='pdf',width=20,height=15)


#lepi.b=tree.fig(group="Lepidopteran",family='ABCBF',node.annot=list(c('PapPol','PluXyl')),label.annot=c('Lepidoptera'))
#ggsave(plot=lepi.b,filename='./CAFE/CAFE_figures/Lepidoptera_ABCBF.pdf',device='pdf',width=20,height=15)


#arac.h=tree.fig(group = 'Arachnid',family = 'ABCH')
#ggsave(plot=arac.h,filename='./CAFE/CAFE_figures/Arac_H.pdf',device='pdf',width=20,height=15)


#arth.h=tree.fig(group = 'Arthropod',family = 'ABCH')
#ggsave(plot=arth.h,filename='./CAFE/CAFE_figures/Arth_H.pdf',device='pdf',width=20,height=15)


