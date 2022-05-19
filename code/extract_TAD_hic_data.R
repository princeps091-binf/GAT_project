library(GenomicRanges)
library(tidyverse)
library(mgcv)
library(readr)
library(furrr)
library(vroom)
options(scipen = 999999999)
options(future.globals.maxSize=20971520000)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
#TAD_file<-"./data/k562_TAD_with_IDs.txt"
#HiC_folder<-"~/Documents/multires_bhicect/data/K562/"

TAD_file<-"~/data_transfer/k562_TAD_with_IDs.txt"
HiC_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/K562/"

#-----------------------------------------
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = 10)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  
#-----------------------------------------
dat_res<-"5kb"
K562_TAD_tbl<-read_delim(TAD_file,col_names = F,delim="\t")
K562_TAD_tbl<-K562_TAD_tbl %>% 
  dplyr::rename(chr=X1,start=X2,end=X3,TAD.ID=X4)

TAD_chr_set<-unique(K562_TAD_tbl$chr)
chr_set<- TAD_chr_set[TAD_chr_set %in% str_split_fixed(list.files(paste0(HiC_folder,dat_res,"/")),pattern = "\\.",2)[,1]] 
for(chromo in chr_set[-1]){
  message("load data for: ", chromo)
  chr_TAD_tbl<-K562_TAD_tbl %>% 
    filter(chr==chromo)
  chr_dat<-compute_chr_res_zscore_fn(HiC_folder,dat_res,chromo,res_num)
  
  tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))
  chr_bin_GRange<-GRanges(seqnames=chromo,
                          ranges = IRanges(start=tmp_bins,
                                           end=tmp_bins + res_num[dat_res] -1
                          ))

  chr_TAD_GRange<-GRanges(seqnames=chr_TAD_tbl$chr,
                          ranges = IRanges(start=chr_TAD_tbl$start,
                                           end=chr_TAD_tbl$end
                          ))
  
  TAD_bin_content_tbl<-findOverlaps(chr_TAD_GRange,chr_bin_GRange) %>% 
    as_tibble %>% 
    mutate(TAD.ID=chr_TAD_tbl$TAD.ID[queryHits],bin=tmp_bins[subjectHits]) %>% 
    group_by(TAD.ID) %>% 
    summarise(bins=list(unique(bin)))
  message("subset TAD data for: ", chromo)
  
  TAD_bin_dat<-chr_dat %>% 
    filter(X1 %in% unique(unlist(TAD_bin_content_tbl$bins)) & X2 %in% unique(unlist(TAD_bin_content_tbl$bins)))
  plan(multisession,workers=15)
  TAD_hic_dat_tbl<-TAD_bin_content_tbl %>% 
#    dplyr::slice_head(n=3) %>% 
    mutate(hic=future_map(bins,function(x){
      TAD_bin_dat %>% 
        filter(X1 %in% x & X2 %in% x)
    }))
  plan(sequential)
  message("saving TAD tables for:", chromo )
  for(i in 1:nrow(TAD_hic_dat_tbl)){
    tmp_tbl<-TAD_hic_dat_tbl$hic[[i]]
    tmp_tbl<-tmp_tbl %>% 
      mutate(bin_a=paste(chr,X1,X1 + res_num[res]-1,sep="_"),
             bin_b=paste(chr,X2,X2 + res_num[res]-1,sep="_")) %>% 
      dplyr::select(bin_a,bin_b,X3,zscore) %>% 
      dplyr::rename(HiC=X3)
    write_tsv(x = tmp_tbl,paste0("/storage/mathelierarea/processed/vipin/group/GAT_data/TAD_data/K562/",TAD_hic_dat_tbl$TAD.ID[i],"_hic_tbl.tsv"))
  }

}
