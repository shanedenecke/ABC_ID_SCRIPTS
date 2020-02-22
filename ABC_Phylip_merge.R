
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
#args[1]='/data2/shane/Transporter_ID/ABC_id/CAFE/Ultrametric_tree/Arachnid_og_sequences'
#print(args[1])
setwd(args[1])

seq.list=list()
names.list=list()
for(i in list.files(full.names = T)[grepl('phy',list.files())]){
 phy=fread(i,sep=' ',header=F) %>% arrange(V1)
 seq.list[[i]]=phy$V2
 names.list[[i]]=phy$V1
}

### are all the names of the list in the correct order?
print(length(unique(names.list))==1)

#### Concatanate all sequences 
concat=c()
for(i in 1:length(seq.list)){
  concat=paste0(concat,seq.list[[i]])  
}

### Add names to sequence file
tax_names=names.list[[1]]
combined.full=cbind(tax_names,concat) %>% data.table()

##Annotate sequence file
len=nchar(concat[1])
headstring=paste0(length(tax_names),len)
colnames(combined.full)=as.character(c(dim(combined.full)[1],len))


fwrite(combined.full,'Full_species.phy',sep=' ')
