compile_locus_results <- function(hess_out_dir) {
  
  message("compiling results")
  
  cond <- read.table(paste0(hess_out_dir,"/COND_Step2"), header=T, stringsAsFactors = F)
  ncond <- read.table(paste0(hess_out_dir,"/NCOND_Step2"), header=T, stringsAsFactors = F)
  
  if(!(all(paste0(cond$chr,cond$start,cond$end) == paste0(ncond$chr,ncond$start,ncond$end)))) {
    stop("cond and noncond file contain different SNPs")
    }
  
  diff <- round(cond$local_h2g / ncond$local_h2g, digits = 3)
  write.table(
    data.frame("chr" = cond$chr, "start" = cond$start, "stop" = cond$end,"cond" = cond$local_h2g, "cond_var" = cond$var, "ncond" = ncond$local_h2g, "ncond_var" = ncond$var, "diff" = diff),
    file = paste0(hess_out_dir,"/compiled_results.csv"),
    quote = F,
    row.names = F,
    sep = ","
  )
}
