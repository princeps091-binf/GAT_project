library(tidyverse)
library(furrr)
library(vroom)

#----------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#----------------------------
fantom_file<-"./data/fantom_ann_entrez_tbl.Rda"
TF_tbl_file<-"./data/TFs_list.tsv"
FANTOM5_tpm_tbl<-vroom("./data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",delim="\t",comment = "#",col_select = 1:7)
TF_tbl<-read_tsv(TF_tbl_file)
FANTOM5_peak_transcript_tbl<-get_obj_in_fn(fantom_file)
K562_CAGE_file<-"./data/K562_CAGE_tbl.Rda"
K562_CAGE_tbl<-get_obj_in_fn(K562_CAGE_file)

FANTOM5_entrez_tbl<-FANTOM5_tpm_tbl %>% 
  dplyr::select(`00Annotation`,entrezgene_id) %>% 
  mutate(entrez.id1=str_split_fixed(entrezgene_id,":",2)[,2]) %>% 
  dplyr::rename(name=`00Annotation`) %>% 
  filter(entrez.id1!="") %>% 
  group_by(name) %>% 
  summarise(entrez.id.l1=list(unique(entrez.id1))) %>% 
  full_join(.,FANTOM5_peak_transcript_tbl  %>% 
              filter(!(is.na(entrez.id))) %>% 
              group_by(name) %>% 
              summarise(entrez.id.l2=list(unique(entrez.id))))
FANTOM5_entrez_tbl<-FANTOM5_tpm_tbl %>% 
  dplyr::select(`00Annotation`,entrezgene_id) %>% 
  mutate(entrez.id1=str_split_fixed(entrezgene_id,":",2)[,2]) %>% 
  dplyr::rename(name=`00Annotation`) %>% 
  filter(entrez.id1!="")

TF_FANTOM5_entrez_tbl<-FANTOM5_entrez_tbl %>% 
  filter(entrez.id1 %in% TF_tbl$entrezgene)

K562_TF_tbl<-K562_CAGE_tbl %>% 
  dplyr::rename(name=`00Annotation`) %>% 
  inner_join(.,TF_FANTOM5_entrez_tbl,by=c("name")) %>% 
  filter(m>0) %>% 
  dplyr::select(name,entrez.id1,m) %>% 
  left_join(.,TF_tbl %>% mutate(entrezgene=as.character(entrezgene)),by=c("entrez.id1"="entrezgene"))
unique(K562_TF_tbl$Symbol)

K562_active_TF_tbl<-K562_TF_tbl %>% 
  distinct(entrez.id1,Symbol) %>% 
  dplyr::rename(entrez.id=entrez.id1)

write_tsv(K562_active_TF_tbl,col_names = T,file = "./data/K562_active_TF.tsv")

K562_TF_tbl %>% 
  ggplot(.,aes(m))+geom_density()+scale_x_log10()
