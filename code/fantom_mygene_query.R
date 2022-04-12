library(tidyverse)
library(furrr)
library(httr)

#----------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#----------------------------
#fantom_file<-"~/data_transfer/FANTOM_annotation_tbl.Rda"
fantom_file<-"./data/FANTOM_annotation_tbl.Rda"

out_file<-"~/data_transfer/fantom_ann_entrez_tbl.Rda"

FANTOM5_peak_transcript_tbl<-get_obj_in_fn(fantom_file)
unique_alt_tbl<- FANTOM5_peak_transcript_tbl %>%
  distinct(alt.name)

chunk_size<-2e3
max_set<-length(unique_alt_tbl$alt.name) - length(unique_alt_tbl$alt.name) %% chunk_size
ends<-seq(chunk_size,max_set,by=chunk_size)
gene_windows<-rbind(cbind(c(1,ends[-length(ends)] +1 ), ends),c(ends[length(ends)]+1,length(unique_alt_tbl$alt.name)))

entrez_set_l<-vector('list',nrow(gene_windows))
for(i in 1:nrow(gene_windows)){
  message(round(i/length(entrez_set_l),digits=4))
  plan(multisession,workers=15)
  tmp_l<-future_map(unique_alt_tbl$alt.name[gene_windows[i,1]:gene_windows[i,2]],function(x){
    tmp<-content(GET(paste0("https://mygene.info/v3/query?q=",x)))$hits
    if(length(tmp)>0){
      return(tmp[[1]]$entrezgene)
    }   else {return(NA_character_)} 
  })
  plan(sequential)
  entrez_set_l[[i]]<-unlist(lapply(tmp_l,function(x)ifelse(is.null(x),NA,x)))
  rm(tmp_l)
}

unique_alt_tbl<-unique_alt_tbl %>% 
  mutate(entrez.id=unlist(entrez_set_l))
FANTOM5_peak_transcript_tbl<-FANTOM5_peak_transcript_tbl %>% 
  left_join(.,unique_alt_tbl)
save(FANTOM5_peak_transcript_tbl,file=out_file)
