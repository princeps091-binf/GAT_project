library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


# Within the /storage/scratch/rafael/GANs/data/familial_BSs/Vipin_by_chr folder containing the TFBS family sites for hg19

bed_file_set<-grep("^chr[0-9]+.bed",list.files(),value=T)

k562_TAD<-import.bed("~/data_transfer/k562_TAD_with_IDs.txt")

tmp_bed<-import.bed(bed_file_set[1])