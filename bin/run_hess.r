
run_hess_step1 <- function(chr, h2g, `reference-panel`, `legend-file`, `partition-file`, out) {
  
  message("run_hess_step1")
  
  command <- paste0(
    "python2.7 ~/Tools/hess-0.3-beta/hess.py",
    " --chrom ", chr,
    " --h2g ", h2g,
    " --reference-panel ", `reference-panel`,
    " --legend-file ", `legend-file`,
    " --partition-file ", `partition-file`,
    " --out ", out
  )
  
  message(command)
  
  system(command)
}


run_hess_step2 <- function(prefix, k=50, out) {
  
  message("run_hess_step2")
  
  command <- paste0(
    "python2.7 ~/Tools/hess-0.3-beta/hess.py",
    " --prefix ", prefix,
    " --k ", k,
    " --out ", out
  )
  
  message(command)
  system(command)
  
}
