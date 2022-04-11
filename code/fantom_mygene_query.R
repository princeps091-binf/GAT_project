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
fantom_file<-"~/data_transfer/FANTOM_annotation_tbl.Rda"
out_file<-"~/data_transfer/fantom_ann_entrez_tbl.Rda"

FANTOM5_peak_transcript_tbl<-get_obj_in_fn(fantom_file)

plan(multisession,workers=15)
tmp_entrez<-future_map_chr(FANTOM5_peak_transcript_tbl$alt.name,function(x){
  tmp<-content(GET(paste0("https://mygene.info/v3/query?q=",x)))$hits
  if(length(tmp)>0){
    return(tmp[[1]]$entrezgene)
  } else{return(NA)}
})
plan(sequential)

FANTOM5_peak_transcript_tbl<-FANTOM5_peak_transcript_tbl %>% 
  mutate(entrez.id=tmp_entrez)

save(FANTOM5_peak_transcript_tbl,file=out_file)
