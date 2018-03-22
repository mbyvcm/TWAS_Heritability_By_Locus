require(argparse, quietly = T, warn.conflicts = F)

source("./bin/append_conditional_summary_stats.r")
source("./bin/run_hess.r")
source("./bin/make_copy_sumstats.r")
source("./bin/remove_hla_snps.r")
source("./bin/compile_locus_results.r")

parser <- ArgumentParser()

parser$add_argument("-n","--array_no")
parser$add_argument("-sumstats_dir", "--sumstats_dir")
parser$add_argument("-cond_dir", "--cond_dir")
parser$add_argument("-dir","--root_dir")
parser$add_argument("-part_dir","--part_dir")
parser$add_argument("-hla","--remove_hla_snps",action="store_true")

args <- parser$parse_args()
array_no   <- args$array_no
root_dir   <- args$root_dir
part_dir <- args$part_dir
hla <- args$remove_hla_snps
gwas_summary_stats_dir <- args$sumstats_dir
outdir <- paste0(root_dir,"output/")
cond_dir <- args$cond_dir

print(array_no)

# get locus details
df <- read.table(paste0(root_dir,"loci.txt"), header=T, stringsAsFactors = F, sep = " ")
locus_chr <- df[array_no,"chr"]
locus_start <- df[array_no,"start"]
locus_stop <- df[array_no,"stop"]
gene  <- df[array_no, "gene"]

# the gwas stats for the specified chromosome (plus chr6 if --hla has been specified)
# will be modified, so need to make a fresh copy into temp folder
tmp_sumstats_dir <- paste0("./tmp/sumstats/",gene,"/")
if(dir.exists(tmp_sumstats_dir)) {
  unlink(paste0(tmp_sumstats_dir,"*"))
    } else { dir.create(tmp_sumstats_dir, recursive = T)}

make_temp_sumstats_file(locus_chr = locus_chr, tmp_sumstats_dir = tmp_sumstats_dir, hla = hla, gwas_summary_stats_dir = gwas_summary_stats_dir)

locus_name <- paste(locus_chr,locus_start,locus_stop, sep = "_")
message(locus_name)

hess_out_dir <- paste0(outdir,"/", locus_name,"/", gene,"/")

if (!(dir.exists(hess_out_dir))) {dir.create(hess_out_dir, recursive = T)}

# generate new conditional / non-conditional files for specified chromosome 
append_conditional_summary_stats(
  chr = locus_chr,
  start = locus_start,
  stop = locus_stop,
  gene = gene,
  gwas_summary_stats_dir = gwas_summary_stats_dir,
  cond_dir = cond_dir
)

tests <- c("COND","NCOND")

for (t in tests) {
  
  file.remove(paste0(tmp_sumstats_dir,"/chr",locus_chr))
  
  modified_file <- paste0("./tmp/modified_sumstats/chr",locus_chr,"_",gene,"_",t)
  dest_file     <- paste0(tmp_sumstats_dir,"/chr",locus_chr)
  
  message(paste0(modified_file," -> ",dest_file))
  
  # move modified file into tmp_sumstats_dir
  if (file.exists(modified_file)) {
    file.copy(from = modified_file, to = dest_file)
  } else {stop("modified file cannot be copied in sumstats dir")}
  
  # fun the first step of hess over each chromosome
  for (chr in seq(22)) {
    message(paste0("chr",chr))
    run_hess_step1(
      chr = chr,
      h2g = paste0(tmp_sumstats_dir,"/chr",chr),
      `reference-panel` = paste0("./resources/refpanel/1kg_phase3_chr",chr,"_eur_5pct_gen.txt.gz"),
      `legend-file` = paste0("./resources/refpanel/1kg_phase3_chr",chr,"_eur_5pct_legend.txt.gz"),
      `partition-file` = paste0(part_dir,"/chr",chr,".bed"),
      out = paste0(hess_out_dir,t)
    )
  }
  
  # run the second step of hess
  run_hess_step2(
    prefix = paste0(hess_out_dir, t),
    k = 50,
    out = paste0(hess_out_dir,t,"_Step2")
  )
}

compile_locus_results(hess_out_dir)
