shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(seqinr))



#setwd('/data2/shane/Transporter_ID/ABC_id/')
sterr=function(x){sd(x)/sqrt(length(x))}

ref.table=read.fasta('./Db_build_temp/Only_ABCs.faa',set.attributes = F,as.string = T)
ref.table.len=data.table(name=names(ref.table),length=sapply(ref.table,nchar))
ref.table.len$fam=gsub('^.+__(.+)__.+$','\\1',ref.table.len$name)

### remove one abnormally short ABCC from tcas
ref.table.len=ref.table.len[name!='TriCas__ABCC__TC014382_TcABCC-5M']

model.mins=ref.table.len %>% group_by(fam) %>% summarize(min=min(length)-sd(length)) %>% data.table()
model.mins=rbind(model.mins,data.table(fam='ABC_Unsorted',min=328))

good.total=list()
short.total=list()
for (i in list.files('./preliminary_ABC/proteomes',full.names = T)){
  fa=read.fasta(i,set.attributes = F,as.string = T)
  test.len=data.table(name=names(fa),length=sapply(fa,nchar))
  test.len$fam=gsub('^.+__(.+)__.+$','\\1',test.len$name)
  
  good.sub=list()
  bad.sub=list()
  for(j in 1:nrow(test.len)){
    sub=test.len[j]
    test.fam=sub$fam
    fam.min=model.mins[fam==test.fam]$min
    
    if(sub$length>fam.min){
      good.sub[[j]]=sub
    }else{
      bad.sub[[j]]=sub
    }
  }
  pass=rbindlist(good.sub)
  fail=rbindlist(bad.sub)
  
  good.total[[i]]=pass
  short.total[[i]]=fail
}

final.pass=rbindlist(good.total)
final.fail=rbindlist(short.total)

final.pass$species=gsub('(^[A-z]+)__AB.+__.+$','\\1',final.pass$name)
final.pass=final.pass[!grepl("Unsorted",species)]
final.pass %>% group_by(species) %>% summarize(count=length(fam)) %>% data.table()


final.fail$species=gsub('(^[A-z]+)__AB.+__.+$','\\1',final.fail$name)
final.fail=final.fail[!grepl("Unsorted",species)]
final.fail %>% group_by(species) %>% summarize(count=length(fam)) %>% data.table()


