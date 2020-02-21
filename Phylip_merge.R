
#rm(list=ls())
## SLC length analysis
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

## already got lengths from expression analysis

args = commandArgs(trailingOnly=TRUE)
#args[1]='/data2/shane/Transporter_ID/ABC_id/CAFE/Ultrametric_tree/Hemipteran_og_sequences'
#print(args[1])
setwd(args[1])

l=list()
l2=list()
for(i in list.files(full.names = T)[grepl('phy',list.files())]){
 phy=fread(i,sep=' ',header=F) 
 l[[i]]=phy$V2
 l2[[i]]=phy$V1
  
  #base=str_extract(i,'[0-9]+')
  #phy=fread(paste(i,'/',base,'_renamed.fasta.aln.trimm.phy',sep=''),sep = ' ',header=F)
  #l[[i]]=phy$V2
  #l2[[i]]=phy$V1
}

full=sapply(l,paste0)

#### No idea what this does 
vec=c()
for(i in 1:length(l)){
  vec=paste0(vec,l[[i]])  
}
tax_names=l2[[i]]
x=cbind(tax_names,vec) %>% data.table()
len=nchar(vec[1])
headstring=paste0(length(tax_names),len)
colnames(x)=as.character(c(dim(x)[1],len))


fwrite(x,'Full_species.phy',sep=' ')
