remove_hla_snps <- function(gwas_summary_stats_dir, tmp_sumstats_dir) {
  
  message("removing HLA region SNPs")
  
  # move chr6
  file.copy(from = paste0(gwas_summary_stats_dir,"/chr",6), to = tmp_sumstats_dir, overwrite = T)
  
  #remove HLA
  chr6_dat <- read.table(paste0(tmp_sumstats_dir,"/chr6"), header = T)
  chr6_dat <- chr6_dat[!(chr6_dat$pos >28e6 & chr6_dat$pos <34e6),]
  write.table(chr6_dat, file =  paste0(tmp_sumstats_dir,"/chr6"), row.names = F, quote = F)
}