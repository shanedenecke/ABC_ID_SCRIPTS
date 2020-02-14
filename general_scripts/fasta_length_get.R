shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(seqinr))
shhh(library(readr))

args = commandArgs(trailingOnly=TRUE)

#args[1]='/data2/shane/Transporter_ID/ABC_id/intermediate/Db_build_temp/Only_ABCs.faa'

fa=read.fasta(args[1],set.attributes = F,as.string = T)

len=sapply(fa,nchar)

final.table=data.table(name=names(fa),length=len)

cat(format_tsv(final.table))
