append_conditional_summary_stats <- function(chr, start, stop, gene, gwas_summary_stats_dir, cond_dir = "./") {
  
  if(!(dir.exists("./tmp/modified_sumstats/"))) {dir.create("./tmp/modified_sumstats/", recursive = T)}
  
chr_stats <- read.table(
  paste0(gwas_summary_stats_dir,"chr",chr),
  header = T,
  stringsAsFactors = F
  ) 
  
locus_stats <- chr_stats[chr_stats$pos > start & chr_stats$pos < stop,]
  
gene_file <- paste0(
  cond_dir,
  "/",
  gene,
  ".loc_1.cond"
  )
  
cond_stats <- read.table(
  gene_file,
  header = T,
  stringsAsFactors = F
  )

i <- match(locus_stats$rsID, cond_stats$SNP)

message(paste0("Number of SNPs in conditional file matched to locus: ", sum(i, na.rm = T)))
    
locus_stats$Z.score <- cond_stats[i,"GWAS_cond.Z"]
locus_stats_cond <- locus_stats[!(is.na(locus_stats$Z.score)),]
    
locus_stats$Z.score <- cond_stats[i, "GWAS.Z"]
locus_stats_noncond <- locus_stats[!(is.na(locus_stats$Z.score)),]
    
    cond_out <- rbind(
      chr_stats[!(chr_stats$pos > start & chr_stats$pos < stop),],
      locus_stats_cond
      )
    
    noncond_out <- rbind(
      chr_stats[!(chr_stats$pos > start & chr_stats$pos < stop),],
      locus_stats_noncond
    )
    
    cond_out <- cond_out[order(cond_out$pos),]
    write.table(cond_out, paste0("./tmp/modified_sumstats/chr",chr,"_",gene,"_COND"), quote = F, row.names = F)
    
    noncond_out <- noncond_out[order(noncond_out$pos),]
    write.table(noncond_out, paste0("./tmp/modified_sumstats/chr",chr,"_",gene,"_NCOND"), quote = F, row.names = F)
}
