library(tidyverse)
library(furrr)
library(vroom)
library(httr)
#
load("./data/K562_CAGE_tbl.Rda")
load("./data/CAGE_tss_ann.Rda")

FANTOM5_tpm_tbl<-vroom("./data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",delim="\t",comment = "#",col_select = 1:7)
#
cage_k562_on<-cage_k562 %>% 
  filter(m>0)

peak_FANTOM_gene_map_tbl<-tot_tss_ann %>% 
  dplyr::select(namess,gene) %>% 
  unnest(cols=c(gene))

FANTOM5_tpm_tbl %>% 
  mutate(entrez.id=str_split_fixed(entrezgene_id,":",2)[,2]) %>% 
  dplyr::select(`00Annotation`,entrez.id) %>% 
  distinct(entrez.id)

FANTOM5_refseq_tbl<-FANTOM5_tpm_tbl %>%
  mutate(transcripts=str_split(association_with_transcript,",")) %>% 
  dplyr::select(`00Annotation`,transcripts) %>% 
  unnest(cols=c(transcripts)) %>% 
  filter(grepl("NM|NR",transcripts))%>% 
  mutate(alt.name=str_extract(transcripts,"NM_[0-9]*|NR_[0-9]*")) %>% 
  dplyr::select(`00Annotation`,alt.name) %>% distinct() %>% 
  dplyr::rename(name=`00Annotation`)



FANTOM5_ENST_tbl<-FANTOM5_tpm_tbl %>%
  mutate(transcripts=str_split(association_with_transcript,",")) %>% 
  dplyr::select(`00Annotation`,transcripts) %>% 
  unnest(cols=c(transcripts)) %>% 
  filter(grepl("ENST",transcripts)) %>% 
  mutate(alt.name=str_extract(transcripts,"ENST[0-9]*")) %>% 
  dplyr::select(`00Annotation`,alt.name) %>% distinct() %>% 
  dplyr::rename(name=`00Annotation`)



FANTOM5_uc_tbl<-FANTOM5_tpm_tbl %>%
  mutate(transcripts=str_split(association_with_transcript,",")) %>% 
  dplyr::select(`00Annotation`,transcripts) %>% 
  unnest(cols=c(transcripts)) %>% 
  filter(grepl("uc",transcripts))%>% 
  mutate(alt.name=str_extract(transcripts,"uc[0-z]+\\.[0-9]+")) %>% 
  dplyr::select(`00Annotation`,alt.name) %>% distinct() %>% 
  dplyr::rename(name=`00Annotation`)


FANTOM5_other_tbl<-FANTOM5_tpm_tbl %>%
  mutate(transcripts=str_split(association_with_transcript,",")) %>% 
  dplyr::select(`00Annotation`,transcripts,association_with_transcript) %>% 
  unnest(cols=c(transcripts)) %>% 
  filter(!(grepl("uc|ENST|NM|NR",transcripts))) %>% 
  filter(!(is.na(transcripts))) %>% 
  mutate(alt.name=str_extract(transcripts,"[A-Z]+[0-9]+")) %>% 
  dplyr::select(`00Annotation`,alt.name) %>% distinct() %>% 
  dplyr::rename(name=`00Annotation`)

FANTOM5_peak_trnscript_tbl<-FANTOM5_other_tbl %>% 
  full_join(.,FANTOM5_uc_tbl) %>% 
  full_join(.,FANTOM5_ENST_tbl) %>% 
  full_join(.,FANTOM5_refseq_tbl) %>% 
  filter(!(is.na(alt.name)))
save(FANTOM5_peak_trnscript_tbl,file = "./data/FANTOM_annotation_tbl.Rda")
plan(multisession,workers=5)
tmp_entrez<-future_map_chr(FANTOM5_peak_trnscript_tbl$alt.name,function(x){
  tmp<-content(GET(paste0("https://mygene.info/v3/query?q=",x)))$hits
  if(length(tmp)>0){
    return(tmp[[1]]$entrezgene)
  } else{return(NA)}
})
plan(sequential)

FANTOM5_peak_trnscript_tbl<-FANTOM5_peak_trnscript_tbl %>% 
  mutate(entrez.id=tmp_entrez)
