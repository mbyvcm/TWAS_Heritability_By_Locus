make_temp_sumstats_file <- function(locus_chr, tmp_sumstats_dir, hla, gwas_summary_stats_dir) {
  
  source("./bin/remove_hla_snps.r")
  
  files <- list.files(gwas_summary_stats_dir, full.names = T)
  # tmp_sumstats_dir <- paste0("./tmp/sumstats/",gene,"/")
  # if(!(dir.exists(tmp_sumstats_dir))){ dir.create(tmp_sumstats_dir, recursive = T)}
  
  file.copy(from = paste0(gwas_summary_stats_dir,"/chr",locus_chr), to = tmp_sumstats_dir, overwrite = T)
  
  if (hla) {
    remove_hla_snps(gwas_summary_stats_dir, tmp_sumstats_dir)
    
    unaltered_files <- files[!(basename(files) %in% c("chr6",paste0("chr",locus_chr)))]
    file.copy(from = unaltered_files, to = tmp_sumstats_dir)
  } else {
    unaltered_files <- files[!(basename(files) %in% paste0("chr",locus_chr))]
    file.copy(from = unaltered_files, to = tmp_sumstats_dir)
    }
}
