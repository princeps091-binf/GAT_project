library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(furrr)
k562_TAD<-import.bed("./data/k562_TAD_with_ID_Family_TFBS_chr1_count.txt")

chr1_TAD<-k562_TAD[seqnames(k562_TAD)=="chr1"]

chr1_TAD %>% as.data.frame %>% 
  mutate(width=abs(end-start)) %>% 
  ggplot(.,aes(width,score))+geom_density_2d_filled(show.legend = F)+scale_x_log10()+scale_y_log10()
Dnase_Grange<-import.bed("./../data/Epi/K562/ENCSR000EKS_rep1_1_rep1_2_rep1_3_se_bwa_biorep_filtered_hotspots_DNase.bed")
chr1_Dnase<-Dnase_Grange[seqnames(Dnase_Grange)=="chr1"]
tmp_opt<-furrr_options(packages = c("GenomicRanges", "IRanges"))


test<-tibble(as.data.frame(findOverlaps(chr1_TAD,chr1_Dnase))) %>% 
  dplyr::rename(query=queryHits,subject=subjectHits) %>% 
  mutate(TAD.id=mcols(chr1_TAD)$name[query])
inter_Grange<-pintersect(chr1_TAD[test$query],chr1_Dnase[test$subject])
export.bed(inter_Grange,con = "./data/chr1_TAD_DNase_open.bed")

k562_TAD_dnase_count_tbl<-read_tsv("./data/k562_TAD_with_ID_Family_TFBS_chr1_dnase_count.txt",col_names = F)

k562_TAD_dnase_count_tbl %>% 
  group_by(X4) %>% 
  summarise(full.count=unique(X5),dnase.count=sum(X7)) %>% 
  ggplot(.,aes(full.count))+
  geom_density(color="red")+
  geom_density(aes(dnase.count),color="blue")+
  scale_x_log10()
ggsave("./img/Dnase_filter_TFBS_count.png")
