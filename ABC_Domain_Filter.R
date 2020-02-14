shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

args = commandArgs(trailingOnly=TRUE)

#setwd('/data2/shane/Transporter_ID/ABC_id')

ipscan=fread('./Filter/IPSCAN.tsv',sep='\t',fill=T)
key=data.table(family=gsub('_','',readLines('./ABC_REF/Input_files/ABC_families.txt')),domains=c(2,2,1,2,1,1,1,1,1))

ipscan$fam=gsub('^.+__(.+)__.+$','\\1',ipscan$V1)







if(nrow(abc.table)>0){
  colnames(abc.table)=c('code','name')
  cat(format_csv(abc.table))
}
