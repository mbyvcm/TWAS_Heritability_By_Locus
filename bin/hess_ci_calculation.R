
outdir <- commandArgs(trailingOnly = T)[1]

cond_files  <- list.files(outdir, pattern = "^COND_Step2$", recursive = T, full.names = T)
ncond_files <- list.files(outdir, pattern = "NCOND_Step2$", recursive = T, full.names = T)

print(cond_files)
print(ncond_files)

calc_95CI <- function(cond_mu, cond_var, main_mu, main_var) {
  
  print(cond_mu)

  RM <- cond_mu
  SR <- cond_var
  SM <- main_mu
  SS <- main_var
  
  if (cond_mu < 0) {cond_mu <- 0}
  
  UM <- 1 - (RM/SM)
  
  VR1 <- (RM^2)/(SM^2)
  VR2 <- (SR/(RM^2)) + (SS/(SM^2))
  VR3 <- (SR/(RM^2))+(SS/(SM^2))-(sqrt(SS*SR)/(SM*RM))
  
  SD=sqrt(VR1*VR2)
  SE=sqrt(VR1*VR3)
  
  UL=round(UM-(1.96*SE), digits = 3)
  UH=round(UM+(1.96*SE), digits = 3)
  
  LL=round(UM-(1.96*SD), digits = 3)
  LH=round(UM+(1.96*SD), digits = 3)
  
  c(SD,SE,UL,UH,LL,LH)
}


out <- lapply(seq(length(cond_files)), function(n) {
  
  coord <- unlist(stringr::str_split(cond_files[n], pattern = "/"))[4]
  chr  <- unlist(stringr::str_split(coord, pattern = "_"))[1]
  start  <- unlist(stringr::str_split(coord, pattern = "_"))[2]
  end  <- unlist(stringr::str_split(coord, pattern = "_"))[3]

  print(coord)

  gene <- unlist(stringr::str_split(cond_files[n], pattern = "/"))[5]
  message(gene)
  
  condDf <- read.table(cond_files[n], header=T, stringsAsFactors = F)
  condDfLocus <- condDf[condDf$chr == chr & condDf$start == start,]
  cond_local_h2g <- condDfLocus$local_h2g
  cond_var       <- condDfLocus$var
  
  ncondDf <- read.table(ncond_files[n], header=T, stringsAsFactors = F)
  ncondDfLocus <- ncondDf[ncondDf$chr == chr & ncondDf$start == start,]
  ncond_local_h2g <- ncondDfLocus$local_h2g
  ncond_var       <- ncondDfLocus$var
  
  diff <- 1 -(cond_local_h2g / ncond_local_h2g)

  ci <- calc_95CI(cond_mu = cond_local_h2g, cond_var = cond_var, main_mu = ncond_local_h2g, main_var = ncond_var)
  
  for (i in c(3,4,5,6)) {
    if (ci[i] < 0) {ci[i] = 0}
    if (ci[i] > 1) {ci[i] = 1}
  }
  
  lower <- min(ci[c(3,5)])
  upper <- max(ci[c(4,6)]) 
  
  out <- data.frame(
    "GENE" = gene,
    "LOCUS_CHR" = chr,
    "LOCUS_START" = start,
    "LOCUS_STOP" = end,
    "MAIN_H2G" = ncond_local_h2g,
    "MAIN_VAR" = ncond_var,
    "COND_H2G" = cond_local_h2g,
    "COND_VAR" = cond_var,
    "DIFF_H2G" = diff,
    "SE_0" = ci[1],
    "SE_1" = ci[2],
    "95CI_L" = lower,
    "95CI_U" = upper
    )
  
  return(out)
  
})

out <- do.call(rbind,out)
out <- out[order(out$LOCUS_CHR, out$LOCUS_START ),]

out <- out[!out$MAIN_H2G < 0,]

out[out$COND_H2G < 0,"COND_H2G"] <- 0

write.table(out, "./hess_output.csv", row.names = F, sep = ",")
