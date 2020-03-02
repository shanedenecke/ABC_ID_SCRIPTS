shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(ggplot2))



############## Directories
setwd('/data2/shane/Transporter_ID/ABC_id')
dir.create('Final_outputs')
dir.create('./Final_outputs/Group_Comparisons')
dir.create('./Final_outputs/Group_Comparisons/ABC_plots')
dir.create('./Final_outputs/Benchmark')
dir.create('./Final_outputs/ABC_fasta')

############## Functions
shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything()) %>% data.table()
  return(b)
}

#### common variables
meta=fread('./ABC_REF/species_metadata/Arthropod_species_metadata.tsv')
counts=fread('./CAFE/ABC_COUNTS_CAFE_FULL.tsv') %>% select(-Desc) %>% rename(fam=`Family ID`)
benchmark.raw=fread('./ABC_REF/species_metadata/ABC_benchmark_counts.csv') %>% select(-Order)



##### Write out summary dictionaries and files


#### Counts
trans=shane.transpose(counts,fam) %>% rename(abbreviation=newcol)
full.counts=merge(meta,trans,by='abbreviation') %>% select(-Common_name) %>% data.table()
full.counts=data.table(select(full.counts,abbreviation:Vory),apply(select(full.counts,ABCA:ABCH),2,as.numeric))
full.counts$ABC_total=rowSums(select(full.counts,ABCA:ABCH))
fwrite(full.counts,'./Final_outputs/Transposed_counts.csv')

#### Copy files and dictionaries from Filter
system('cp -r /data2/shane/Transporter_ID/ABC_id/Filter/Full_transporters/* /data2/shane/Transporter_ID/ABC_id/Final_outputs/')
system('cat /data2/shane/Transporter_ID/ABC_id/Final_outputs/proteomes/* > /data2/shane/Transporter_ID/ABC_id/Final_outputs/proteomes/Combined_ABC_proteomes.faa')
dict.list=list()
for(i in list.files('./Filter/Full_transporters/dicts',full.names = T)){
  dict.list[[i]]=fread(i)
}
fwrite(rbindlist(dict.list),'./Final_outputs/dicts/Combined_total_dictionary.tsv',sep='\t')
#### combine dictionaries



################# Benchmarking against known datasets 
merged.benchmark=select(full.counts,Species_name,ABC_total) %>% merge(benchmark.raw,by='Species_name') %>% rename(ABC_This_Study=ABC_total)
merged.benchmark$Difference=merged.benchmark$ABC_This_Study-merged.benchmark$Lit_ABC_count
merged.benchmark$Percent_Difference=merged.benchmark$Difference/merged.benchmark$ABC_This_Study


gp=ggplot(merged.benchmark,aes(x=Species_name,y=Percent_Difference))
gp=gp+geom_bar(stat='identity')
gp=gp+geom_hline(yintercept=0,linetype=1,size=2,color='red')
gp=gp+scale_y_continuous(breaks=seq(-.5,.5,by=.1),limits=c(-.5,.5))
gp=gp+theme_bw()
gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
            axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
            strip.text=element_text(size=20),strip.background=element_rect("white"),
            axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1),
            legend.position = 'none',plot.title = element_text(hjust = 0.5))

print(gp)
ggsave(plot=gp,filename='./Final_outputs/Benchmark/Benchmark_graph.pdf')



###################### Perform ANOVA on groups and looks at plots
iter.i=select(full.counts,ABCA:ABCH) %>% colnames()
iter.j=iter.j=select(full.counts,Taxonomic_Classification:Vory) %>% colnames()

anova.l=list()
for(i in iter.i){
  for(j in iter.j){
    
    #### Filter Data for those groups which have over 5 individuals
    good=names(table(full.counts[[j]]))[table(full.counts[[j]])>5]
    sub=full.counts[full.counts[[j]] %in% good]
    
    
    ### Make ggplot
    gp=ggplot(data=sub,aes_string(x=j,y=i,fill=j))
    gp=gp+geom_boxplot(outlier.size=1)
    gp=gp+theme_bw()
    gp=gp+theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),
                axis.ticks.x=element_line(),panel.border=element_rect(colour="black",fill=NA),
                strip.text=element_text(size=20),strip.background=element_rect("white"),
                axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1),
                legend.position = 'none',plot.title = element_text(hjust = 0.5))
    
    #print(gp)
    ggsave(filename=paste0('./Final_outputs/Group_Comparisons/ABC_plots/',i,'_',j,'.png'),gp,device='png',height=10,width=10)
    
    
    ##### ANOVA
    model=aov(formula=sub[[i]]~sub[[j]])
    pval=summary(model)[[1]][["Pr(>F)"]][1]
    ef=max(abs(model$coefficients[!grepl('Intercept',names(model$coefficients))]))
    
    tuk=TukeyHSD(model)[[1]]
    comp=rownames(tuk)
    tuk=data.table(tuk) %>% rename(pval=`p adj`)
    tuk$comparison=comp
    
    anova.l[[paste(i,j,sep='_')]]=data.table(tuk,family=i,co_variable=j) 
  }
}



anova.summary=rbindlist(anova.l)  
anova.summary$bonf=p.adjust(anova.summary$pval,method='bonferroni') 
anova.filter=anova.summary %>% arrange(bonf) %>% filter(bonf<1e-2)  %>% data.table()
anova.distinct=anova.filter %>% select(family,co_variable) %>% unique.data.frame()
fwrite(anova.filter,'./Final_outputs/Group_Comparisons/ABC_plots')

