shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(readr))
shhh(library(stringr))

args = commandArgs(trailingOnly=TRUE)

#setwd('/data2/shane/Transporter_ID/ABC_id')
#args[1]='./ABC_search/MyzPer/total_ABC_recip_blast.tsv'


abc.number=function(x){
  l=list()
  for(i in unique(x$family)){
    sub=x[family==i]
    sub$family=paste0(sub$family,'_',1:nrow(sub))
    l[[i]]=sub
  }
  return(rbindlist(l)$family)
}



recip_blast=fread(args[1]) %>% select(V1,V2,V4) %>% rename(query=V1,subject=V2,evalue=V4)
recip_blast$subject=gsub('^.+__(.+)__.+$','\\1',recip_blast$subject)

species=unlist(strsplit(args[1],'\\/'))[3]
target.fam=as.character(str_match(unlist(strsplit(args[1],'\\/'))[4],'ABC[A-Z]'))
abc.total=list()
filter.list=list()


for(i in unique(recip_blast$query)){
  #i="LOC111041090"
  sub=recip_blast[query==i]
  top_hit_fam=sub[query==i]$subject[1]
  fams=table(sub$subject)
  evalues=sub$evalue
  
  
  if(!grepl('ABC[A-Z]',sub$subject[1])){ ## filter out any that have top blast hit not an SLC
    filter.list[[i]]=data.table(geneid=i,family=top_hit_fam)
  } else if(length(fams[fams>3])>0){ ### Sort into a family if there are 4 or more of the top 5 blast hits in that family
    abc.total[[i]]=data.table(geneid=i,family=names(fams[fams>3]))
  } else if((evalues[2]/evalues[1]> 1e30) & grepl('ABC',top_hit_fam)){ ## Keep where top hit is SLC and overwhelmingly significant
    abc.total[[i]]=data.table(geneid=i,family=top_hit_fam)
  } else { ## keep where 4 out of 5 are SLCs but can be any family
    abc.total[[i]]=data.table(geneid=i,family='ABC_Unsorted_')
  }
}




filter.table=rbindlist(filter.list) 
abc.table.reduce=rbindlist(abc.total) %>% arrange(family) %>% data.table()
abc.table.reduce$family=paste0(species,'__',abc.table.reduce$family,'__',abc.number(abc.table.reduce))



if(nrow(filter.table)>0){
  colnames(filter.table)=c('code','name') 
  #filter.output=filter.table %>% merge(source.table,by='code',all=T)
  #write.csv(filter.output,'./ABC_filtered_out.csv')
}

if(nrow(abc.table.reduce)>0){
  colnames(abc.table.reduce)=c('code','name')
  cat(format_csv(abc.table.reduce))
}
