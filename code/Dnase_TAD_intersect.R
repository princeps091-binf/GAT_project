library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(furrr)
k562_TAD<-import.bed("./data/k562_TAD_with_ID_Family_TFBS_chr1_count.txt")

chr1_TAD<-k562_TAD[seqnames(k562_TAD)=="chr1"]

hist(log10(mcols(chr1_TAD)$score))

Dnase_Grange<-import.bed("./../data/Epi/K562/ENCSR000EKS_rep1_1_rep1_2_rep1_3_se_bwa_biorep_filtered_hotspots_DNase.bed")
chr1_Dnase<-Dnase_Grange[seqnames(Dnase_Grange)=="chr1"]
tmp_opt<-furrr_options(packages = c("GenomicRanges", "IRanges"))

plan(multisession, workers=5)

test<-tibble(as.data.frame(findOverlaps(chr1_TAD,Dnase_Grange))) %>% 
  dplyr::rename(query=queryHits,subject=subjectHits) %>% 
#  slice(1:5) %>% 
    mutate(inter.Grange=future_pmap(list(query,subject),.options = tmp_opt,function(query,subject){
    
    return(IRanges::intersect(chr1_TAD[query],chr1_Dnase[subject]))
    
  })) %>% 
  mutate(TAD.id=mcols(chr1_TAD)$name[query])
plan(sequential)

inter_Grange<-do.call("c",test$inter.Grange)

mcols(inter_Grange)<-tibble(ID=test$TAD.id)