shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(ape))
shhh(library(ggtree))
shhh(library(ggplot2))
shhh(library(treeio))
library(RColorBrewer)

## Import arguments
args = commandArgs(trailingOnly=TRUE)
#args[1]= '~/Transporter_ID/ABC_id/phylo/clean_trees'
setwd(args[1])

if(length(args)==2){
  species=readLines(args[2])
}else{
    species=''
  }

calc.xmax=function(x){
  x=tbl 
  tips=tbl[grepl('[A-Z]',label)]$node
  bl=c()
  for(i in tips){
    anc=ancestor(base.tree,i)
    total=x[node %in% anc]$branch.length %>% sum(na.rm=T)
    bl=c(bl,total)
  }
  return(max(bl))
}



di2multi4node <- function (phy, tol = 0.5) {
  library(ape)
  # Adapted di2multi function from the ape package to plot polytomies
  # based on numeric node support values
  # (di2multi does this based on edge lengths)
  # Needs adjustment for unrooted trees as currently skips the first edge
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  if (is.na(as.numeric(phy$node.label[2])))
    stop("node labels can't be converted to numeric values")
  if (is.null(phy$node.label))
    stop("the tree has no node labels")
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) 
    return(phy)
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) 
      next
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}






### create iterable list for trees
iter.trees=list.files()[grepl('RAxML_bipartitions.',list.files(),fixed=T)]

for(i in iter.trees){
  b=gsub('RAxML_bipartitions.(.+).[a-z]{3}$','\\1',i)
  base.tree=read.tree(i)
  
  out.g=base.tree$tip.label[grepl('OUT',base.tree$tip.label)]
  
  #if(grepl('OUT',out.g)){
  #  base.tree=root(base.tree,outgroup=out.g,edgelabel = T)
  #}
  
  base.tree$node.label=as.numeric(base.tree$node.label)
  base.tree$node.label[is.na(base.tree$node.label)]=99.99
  col.tree=di2multi4node(base.tree,20)
  col.tree$node.label[col.tree$node.label==99.99]=NA
  
  tbl=as_tibble(col.tree) %>% data.table()
  ## set node colors based on bootstrap values
  node.cols=c()
  for(j in as.numeric(col.tree$node.label)){
    if(is.na(j)){node.cols=c(node.cols,'black')
    }else if(j>80){node.cols=c(node.cols,'greenyellow')
    }else if(j>50){node.cols=c(node.cols,'skyblue')
    }else{node.cols=c(node.cols,'pink')}
  }
  
  ## set tip colors based on species
  unispec=readLines('~/Transporter_ID/ABC_id/ABC_REF/Input_files/Phylo_list.txt')
  unispec=c(unispec,out.g)
  names(unispec)=brewer.pal(length(unispec),'Dark2')
  tip.cols=c()
  for(j in col.tree$tip.label){tip.cols=c(tip.cols,names(unispec)[sapply(unispec,grepl,j)])}
  
  #xmax=calc.xmax(tbl) * 1.5
  
  ### make plot 
  gp=ggtree(col.tree,size=2)
  gp=gp+geom_tiplab(size=6,fontface='bold',col=tip.cols)
  gp=gp+geom_nodepoint(size=8,col=node.cols)
  
  ### set scales
  xmax=layer_scales(gp)$x$range$range[2]*1.3
  h=max(15,length(col.tree$tip.label)/3)
  #w=xmax*7
  w=h*1.5
  if(length(tip.cols)<40){
    hj=w/20
  }else{
    hj=w/30
  }
  
  gp=gp+geom_nodelab(hjust=hj,vjust=.3,size=3,fontface='bold',col='black')
  gp=gp+ggtitle(b)
  gp=gp+theme_tree()
  
  gp=gp+scale_x_continuous(limits=c(0,xmax))
  #gp=gp+scale_x_continuous(limits=c(0,xmax))
  #print(gp)
  ggsave(paste0('./ggrtree_phylogeny',b,'.pdf'),plot=gp,device='pdf',height=h,width=w,limitsize = F)
}