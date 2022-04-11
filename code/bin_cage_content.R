library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(furrr)

chr_dat<-read_tsv("./data/k562_chr1_5kb_KR_.txt",col_names = F)
#
load("./data/K562_CAGE_tbl.Rda")
#
cage_k562_on<-cage_k562 %>% 
  filter(m>0)

tmp_coord<-str_split_fixed(cage_k562_on$`00Annotation`[-(1:2)],pattern = ":|\\.\\.|,",n = 4)

cage_Grange<-GRanges(seqnames=tmp_coord[,1],
                         ranges = IRanges(start=as.numeric(tmp_coord[,2]),
                                          end=as.numeric(tmp_coord[,3])
                         ))

chr1_bin<-unique(c(chr_dat$X1,chr_dat$X2))

bin_Grange<-GRanges(seqnames="chr1",
                    ranges = IRanges(start=as.numeric(chr1_bin),
                                     end=as.numeric(chr1_bin) + 5000 -1
                    ))

tibble(chr="chr1",bin=chr1_bin,cage.count=countOverlaps(bin_Grange,cage_Grange)) %>% 
  mutate(cage=ifelse(cage.count>0,"cage","no cage")) %>% 
  group_by(cage) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(x="bin",n,fill=cage))+geom_bar(stat="identity")
#0.06876055
ggsave("./img/CAGE_bin_count_chr1.png")

k562_TAD<-import.bed("./data/k562_TAD_with_ID_Family_TFBS_chr1_count.txt")

chr1_TAD<-k562_TAD[seqnames(k562_TAD)=="chr1"]

as.data.frame(chr1_TAD) %>% 
mutate(bin.count=countOverlaps(chr1_TAD,bin_Grange),cage.count=countOverlaps(chr1_TAD,cage_Grange))

tmp_coord<-str_split_fixed(cage_k562$`00Annotation`[-(1:2)],pattern = ":|\\.\\.|,",n = 4)

cage_all_Grange<-GRanges(seqnames=tmp_coord[,1],
                     ranges = IRanges(start=as.numeric(tmp_coord[,2]),
                                      end=as.numeric(tmp_coord[,3])
                     ))
tibble(chr="chr1",bin=chr1_bin,cage.count=countOverlaps(bin_Grange,cage_all_Grange)) %>% 
  mutate(cage=ifelse(cage.count>0,"cage","no cage")) %>% 
  group_by(cage) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(x="bin",n,fill=cage))+geom_bar(stat="identity")
#0.1742072
ggsave("./img/CAGE_all_bin_count_chr1.png")
